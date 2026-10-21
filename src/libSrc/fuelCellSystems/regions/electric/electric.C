/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "electric.H"
#include "IOdictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

#include "sigmaModelList.H"
#include "dissolvedModel.H"
#include "activationOverpotentialModel.H"

#include "fuelCellSystem.H"
#include "lduPrimitiveMesh.H"
#include "lduMatrix.H"
#include "lduPrimitiveProcessorInterface.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{
    defineTypeNameAndDebug(electric, 0);

    addToRunTimeSelectionTable
    (
        regionType,
        electric,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionTypes::electric::electric
(
    const fvMesh& mesh,
    const word& regionName
)
:
    regionType(mesh, regionName),

    mesh_(mesh),
    i_
    (
        IOobject
        (
            "i",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedVector(dimCurrent/dimArea, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    phi_
    (
        IOobject
        (
            activationOverpotentialModel::phiName,
            this->time().timeName(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        *this
    ),
    j_
    (
        IOobject
        (
            activationOverpotentialModel::jName,
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimCurrent/dimVol, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaField_
    (
        IOobject
        (
            "sigma",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(j_.dimensions()/phi_.dimensions()*dimArea, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    T_
    (
        IOobject
        (
            "T",
            this->time().timeName(),
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar(dimTemperature, dict_.getOrDefault<scalar>("T0", 293.15)),
        zeroGradientFvPatchScalarField::typeName
    ),
    relax_(dict_.lookupOrDefault<scalar>("relax", 0.0)),
    control_(dict_.lookupOrDefault<Switch>("control", false)),
    dissolveOnOff_(dict_.lookupOrDefault<Switch>("dissolveOnOff", false)),
    patchName_(word::null),
    zoneName_(word::null),
    cellZoneIDs_(),
    galvanostatic_(true),
    ibar_(nullptr),
    voltage_(nullptr),
    coupled_(false),
    solvedByCoupler_(false),
    electronRegion_(word::null),
    omega_(1.0)
{
    if (dissolveOnOff_)
    {
        dissolved_ = dissolvedModel::New(*this, dict_.subDict("dissolved"));
    }

    sigma_.reset
    (
        new sigmaModelList
        (
            *this,
            dict_.subDict("sigma")
        )
    );

    if (control_)
    {
        const dictionary& controlDict = dict_.subDict("galvanostatic");
        patchName_ = controlDict.get<word>("patchName");
        galvanostatic_ = controlDict.get<Switch>("active");

        if (galvanostatic_)
        {
            ibar_.reset(Function1<scalar>::New("ibar", controlDict, this));
        }
        else
        {
            voltage_.reset(Function1<scalar>::New("voltage", controlDict, this));
        }

        if (controlDict.found("zoneName"))
        {
            zoneName_ = controlDict.get<word>("zoneName");
            cellZoneIDs_ = this->cellZones().indices(zoneName_);
        }
    }

    //- Coupled phiE-phiI solve, set in the ion region:
    //  coupled { active yes; electronRegion phiE; omega 1; }
    if (dict_.found("coupled"))
    {
        const dictionary& coupledDict = dict_.subDict("coupled");
        coupled_ = coupledDict.get<Switch>("active");

        if (coupled_)
        {
            electronRegion_ = coupledDict.get<word>("electronRegion");
            omega_ = coupledDict.getOrDefault<scalar>("omega", 1.0);

            //- Their existence switches ButlerVolmer to coupled mode
            dJdEta_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "dJdEta",
                        this->time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    *this,
                    dimensionedScalar(j_.dimensions()/phi_.dimensions(), Zero)
                )
            );
            etaE_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "etaE",
                        this->time().timeName(),
                        *this,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    *this,
                    dimensionedScalar(phi_.dimensions(), Zero)
                )
            );

            //- The electron region (constructed before) is solved from here
            electric& electron = const_cast<electric&>
            (
                refCast<const electric>
                (
                    this->time().lookupObject<regionType>(electronRegion_)
                )
            );
            electron.solvedByCoupler_ = true;

            Info<< "Region " << name() << ": coupled solve with "
                << electronRegion_ << ", omega = " << omega_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionTypes::electric::~electric()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionTypes::electric::solve()
{
    if (solvedByCoupler_)
    {
        return;
    }

    if (coupled_)
    {
        solveCoupled();
        return;
    }

    Info << "\nSolve for region " << name() << ":\n" << endl;

    tmp<fvScalarMatrix> phiEqn
    (
      - fvm::laplacian(sigmaField_, phi_, "laplacian(sigma,phi)")
      - j_
    );

    //- Set reference values
    if (phi_.needReference())
    {
        //- Update the potential field.
        //- Electro-Neutral, the total electron+proton flux should be zero
        const scalarField& source = j_;
        const scalarField& volume = this->V();
        scalarField sum = source*volume;
        scalar iDot(Foam::gSum(sum)/Foam::gSum(volume));

        for (label id = 0; id < nCells(); id++)
        {
            //- The option has to be true (forceReference)
            //- otherwise, errors happen in parallel simulations
            phiEqn->setReference(id, phi_[id] - iDot*relax_, true);
        }
    }

    //phiEqn->relax();

    phiEqn->solve();

    i_ = -sigmaField_ * fvc::grad(phi_);
    i_.correctBoundaryConditions();

    if (dissolved_.valid())
    {
        dissolved_->solve();
    }
}


void Foam::regionTypes::electric::solveCoupled()
{
    //- Called in the ion region. Unknowns [phiE; phiI] in one matrix.
    //- J of both regions holds the base-point current J* (= realised
    //- current of the previous iteration), ButlerVolmer has written
    //- dJdEta = G'(eta*) and etaE = eta* + E in the CL cells.
    electric& eR = const_cast<electric&>
    (
        refCast<const electric>
        (
            this->time().lookupObject<regionType>(electronRegion_)
        )
    );
    electric& iR = *this;

    Info<< "\nCoupled solve for regions " << eR.name() << " + " << iR.name()
        << ":\n" << endl;

    //- Segregated matrices: -laplacian(sigma, phi) - J
    fvScalarMatrix eEqn
    (
      - fvm::laplacian(eR.sigmaField_, eR.phi_, "laplacian(sigma,phi)")
      - eR.j_
    );
    fvScalarMatrix iEqn
    (
      - fvm::laplacian(iR.sigmaField_, iR.phi_, "laplacian(sigma,phi)")
      - iR.j_
    );

    const label nE = eR.nCells();
    const label nI = iR.nCells();
    const label nCells = nE + nI;

    //- Cell pairs: ion cell -> electron cell (via the parent mesh).
    //- In parallel all copies of a cell are on the same processor
    //- (decomposeParID), so the pairs are always local.
    DynamicList<label> pairE;
    DynamicList<label> pairI;
    forAll(iR.cellMapIO_, cellI)
    {
        const auto iter = eR.cellMap_.cfind(iR.cellMapIO_[cellI]);
        if (iter.good())
        {
            pairI.append(cellI);
            pairE.append(iter.val());
        }
    }
    const label nP = pairI.size();

    //- Addressing: phiE faces, phiI faces (offset nE), one link per pair
    const labelUList& lE = eR.lduAddr().lowerAddr();
    const labelUList& uE = eR.lduAddr().upperAddr();
    const labelUList& lI = iR.lduAddr().lowerAddr();
    const labelUList& uI = iR.lduAddr().upperAddr();
    const label nFE = lE.size();
    const label nFI = lI.size();
    const label nF = nFE + nFI + nP;

    labelList lower(nF);
    labelList upper(nF);
    scalarField upperCoeffs(nF);

    const scalarField& upperE = eEqn.upper();
    const scalarField& upperI = iEqn.upper();

    for (label f = 0; f < nFE; ++f)
    {
        lower[f] = lE[f];
        upper[f] = uE[f];
        upperCoeffs[f] = upperE[f];
    }
    for (label f = 0; f < nFI; ++f)
    {
        lower[nFE + f] = nE + lI[f];
        upper[nFE + f] = nE + uI[f];
        upperCoeffs[nFE + f] = upperI[f];
    }

    //- Diagonal and source, boundary contributions folded in
    //- (as fvMatrix::solve does with addBoundaryDiag/addBoundarySource).
    //- Processor patches become interfaces of the combined matrix: their
    //- faceCells are offset like the cells, the communication is done by
    //- the processor patch fields of phi, which take the faceCells from
    //- the combined addressing (lduAddr.patchAddr).
    scalarField diag(nCells);
    scalarField source(nCells);
    scalarField psi(nCells);

    SubList<scalar>(diag, nE, 0) = eEqn.diag();
    SubList<scalar>(diag, nI, nE) = iEqn.diag();
    SubList<scalar>(source, nE, 0) = eEqn.source();
    SubList<scalar>(source, nI, nE) = iEqn.source();
    SubList<scalar>(psi, nE, 0) = eR.phi_.primitiveField();
    SubList<scalar>(psi, nI, nE) = iR.phi_.primitiveField();

    PtrList<const lduInterface> procInterfaces;
    lduInterfaceFieldPtrsList interfaces;
    FieldField<Field, scalar> interfaceBouCoeffs;
    FieldField<Field, scalar> interfaceIntCoeffs;

    auto addBoundary = [&](const electric& r, const fvScalarMatrix& eqn, label off)
    {
        const lduInterfaceFieldPtrsList phiInterfaces =
            r.phi_.boundaryField().scalarInterfaces();

        forAll(r.phi_.boundaryField(), patchi)
        {
            const labelUList& fc = r.boundary()[patchi].faceCells();
            const scalarField& ic = eqn.internalCoeffs()[patchi];
            const scalarField& bc = eqn.boundaryCoeffs()[patchi];

            forAll(fc, k)
            {
                diag[off + fc[k]] += ic[k];
            }

            if (!r.phi_.boundaryField()[patchi].coupled())
            {
                forAll(fc, k)
                {
                    source[off + fc[k]] += bc[k];
                }
                continue;
            }

            if (!isA<processorFvPatch>(r.boundary()[patchi]))
            {
                FatalErrorInFunction
                    << "Coupled patch " << r.boundary()[patchi].name()
                    << " in region " << r.name() << " is not supported"
                    << exit(FatalError);
            }

            const processorFvPatch& pp =
                refCast<const processorFvPatch>(r.boundary()[patchi]);

            labelList fcOff(fc.size());
            forAll(fc, k)
            {
                fcOff[k] = off + fc[k];
            }

            const label k = procInterfaces.size();

            procInterfaces.resize(k + 1);
            procInterfaces.set
            (
                k,
                new lduPrimitiveProcessorInterface
                (
                    fcOff,
                    pp.myProcNo(),
                    pp.neighbProcNo(),
                    pp.forwardT(),
                    pp.tag(),
                    pp.comm()
                )
            );

            interfaces.resize(k + 1);
            interfaces.set(k, phiInterfaces.get(patchi));

            interfaceBouCoeffs.resize(k + 1);
            interfaceBouCoeffs.set(k, new scalarField(bc));
            interfaceIntCoeffs.resize(k + 1);
            interfaceIntCoeffs.set(k, new scalarField(ic));
        }
    };
    addBoundary(eR, eEqn, 0);
    addBoundary(iR, iEqn, nE);

    //- Coupling terms, r = omega*V*G'(eta*):
    //  row e: (aP + r) phiE - r phiI = -V J* + r (eta* + E)
    //  row i: (aP + r) phiI - r phiE = +V J* - r (eta* + E)
    //  (the V J* parts are already in the source through "- J")
    const scalarField& V = iR.V();
    const scalarField& dJdEta = dJdEta_();
    const scalarField& etaE = etaE_();
    scalarField r(nP);

    forAll(pairI, p)
    {
        const label e = pairE[p];
        const label i = pairI[p];
        r[p] = omega_*V[i]*dJdEta[i];

        lower[nFE + nFI + p] = e;
        upper[nFE + nFI + p] = nE + i;
        upperCoeffs[nFE + nFI + p] = -r[p];

        diag[e] += r[p];
        diag[nE + i] += r[p];
        source[e] += r[p]*etaE[i];
        source[nE + i] -= r[p]*etaE[i];
    }

    //- Upper-triangular face order
    const labelList oldToNew
    (
        lduPrimitiveMesh::upperTriOrder(nCells, lower, upper)
    );

    labelList lowerNew(nF);
    labelList upperNew(nF);
    scalarField upperCoeffsNew(nF);
    forAll(oldToNew, f)
    {
        lowerNew[oldToNew[f]] = lower[f];
        upperNew[oldToNew[f]] = upper[f];
        upperCoeffsNew[oldToNew[f]] = upperCoeffs[f];
    }

    lduInterfacePtrsList meshInterfaces(procInterfaces.size());
    forAll(procInterfaces, k)
    {
        meshInterfaces.set(k, &procInterfaces[k]);
    }
    const lduSchedule schedule
    (
        lduPrimitiveMesh::nonBlockingSchedule<processorLduInterface>
        (
            meshInterfaces
        )
    );

    //- The mesh takes ownership of procInterfaces
    lduPrimitiveMesh lduMesh
    (
        nCells,
        lowerNew,
        upperNew,
        procInterfaces,
        schedule,
        UPstream::worldComm
    );

    //- Symmetric matrix: only diag and upper are set
    lduMatrix M(lduMesh);
    M.diag() = diag;
    M.upper() = upperCoeffsNew;

    solverPerformance perf = lduMatrix::solver::New
    (
        "phiCoupled",
        M,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        iR.solverDict("phiCoupled")
    )->solve(psi, source);

    Info<< perf.solverName() << ":  Solving for phiCoupled"
        << ", Initial residual = " << perf.initialResidual()
        << ", Final residual = " << perf.finalResidual()
        << ", No Iterations " << perf.nIterations() << endl;

    //- Scatter back
    eR.phi_.primitiveFieldRef() = SubList<scalar>(psi, nE, 0);
    iR.phi_.primitiveFieldRef() = SubList<scalar>(psi, nI, nE);
    eR.phi_.correctBoundaryConditions();
    iR.phi_.correctBoundaryConditions();

    //- Realised current J = J* + (r/V)(etaNew - eta*), it becomes the
    //- base-point current of the next iteration
    scalar maxDelta = 0;
    scalar Ipos = 0;
    scalar Ineg = 0;

    forAll(pairI, p)
    {
        const label e = pairE[p];
        const label i = pairI[p];

        const scalar delta = eR.phi_[e] - iR.phi_[i] - etaE[i];
        const scalar Jreal = iR.j_[i] + r[p]/V[i]*delta;

        iR.j_[i] = Jreal;
        eR.j_[e] = -Jreal;

        maxDelta = max(maxDelta, mag(delta));
        Ipos += max(Jreal, scalar(0))*V[i];
        Ineg += min(Jreal, scalar(0))*V[i];
    }
    eR.j_.correctBoundaryConditions();
    iR.j_.correctBoundaryConditions();

    label nPairs = nP;
    reduce(nPairs, sumOp<label>());
    reduce(maxDelta, maxOp<scalar>());
    reduce(Ipos, sumOp<scalar>());
    reduce(Ineg, sumOp<scalar>());

    Info<< "Coupled: pairs = " << nPairs
        << ", max|eta_new - eta*| = " << maxDelta
        << " V, I+ (A) = " << Ipos
        << ", I- (A) = " << Ineg
        << ", rel. imbalance = " << (Ipos + Ineg)/max(Ipos, VSMALL)
        << endl;

    //- Current densities and dissolved model, as in the segregated solve
    eR.i_ = -eR.sigmaField_*fvc::grad(eR.phi_);
    eR.i_.correctBoundaryConditions();
    iR.i_ = -iR.sigmaField_*fvc::grad(iR.phi_);
    iR.i_.correctBoundaryConditions();

    if (eR.dissolved_.valid())
    {
        eR.dissolved_->solve();
    }
    if (iR.dissolved_.valid())
    {
        iR.dissolved_->solve();
    }
}


void Foam::regionTypes::electric::setRDeltaT()
{
    //- Do nothing, add if necessary.
}


void Foam::regionTypes::electric::correct()
{
    sigma_->correct(sigmaField_);
    sigmaField_.correctBoundaryConditions();

    if (control_)
    {
        //- Update the potential field.
        const scalarField& source = j_;
        const scalarField& volume = this->V();
        scalar iDot(0.0);

        if (cellZoneIDs_.size())
        {
            forAll(cellZoneIDs_, zoneI)
            {
                const labelList& cells = this->cellZones()[cellZoneIDs_[zoneI]];

                forAll(cells, i)
                {
                    const label cellI = cells[i];

                    iDot += source[cellI] * volume[cellI];
                }
            }

            reduce(iDot, sumOp<scalar>());
        }
        else
        {
            iDot = Foam::gSum(source*volume);
        }

        {
            label patchID = this->boundaryMesh().findPatchID(patchName_);

            if (patchID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find " << patchName_
                    << "please check the patchName" << exit(FatalError);
            }

            scalarField phiBoundary = phi_.boundaryFieldRef()[patchID];

            scalar ibar0 = iDot/Foam::gSum(this->magSf().boundaryField()[patchID]);

            if (galvanostatic_)
            {
                phi_.boundaryFieldRef()[patchID] ==
                    phiBoundary + relax_*(ibar0 - ibar_->value(time().value()));
            }
            else
            {
                phi_.boundaryFieldRef()[patchID] == voltage_->value(time().value());
            }
            
            Info << "ibar: " << ibar0 << "\t"<< "voltage: " << Foam::gAverage(phiBoundary) << endl;
        }
    }

    if (dissolved_.valid())
    {
        dissolved_->correct();
    }
}


void Foam::regionTypes::electric::mapToCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " to Cell " << nl << endl;

    //- heat source
    volScalarField heatSource = (i_ & i_)/sigmaField_;
    volScalarField heatSource0
    (
        IOobject
        (
            "heatSource",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    heatSource0.rmap(heatSource, cellMapIO_);

    fuelCell.Qdot() += heatSource0;
}


void Foam::regionTypes::electric::mapFromCell
(
    fuelCellSystem& fuelCell
)
{
    Info << "Map " << name() << " from Cell " << nl << endl;

    scalarField& T = T_;

    forAll(T, cellI)
    {
        T[cellI] = fuelCell.T()[cellMapIO_[cellI]];
    }

    T_.correctBoundaryConditions();
}

// ************************************************************************* //
