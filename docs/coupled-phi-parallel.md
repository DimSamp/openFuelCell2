# Coupled φ_E–φ_I solve in parallel

Branch `coupled-phi-parallel`, on top of the serial version (`coupled-phi`,
commit `bbc7653`). The only file changed is
`src/libSrc/fuelCellSystems/regions/electric/electric.C`, in
`electric::solveCoupled()`. The method, the case switches and the
ButlerVolmer part are unchanged.

## How to switch between versions

```bash
cd ~/OpenFOAM/cfdlab-v2312/applications/openFuelCell2
git checkout coupled-phi            # serial coupled solve
git checkout coupled-phi-parallel   # serial + parallel coupled solve
git checkout master                 # original openFuelCell2 (segregated only)

# after every checkout, rebuild the library
cd src/libSrc/fuelCellSystems && openfoam2312 -c 'wmake libso'
```

## What changed and why

In serial, the combined matrix holds the phiE cells (0 … N_E−1), then the phiI
cells (N_E … N_E+N_I−1), plus one link per CL cell pair. It had no
interfaces, so any processor patch stopped the run with a FatalError.

Parallel needs only the processor patches of both regions as interfaces of
the combined matrix:

1. **Serial-only check removed.** The `Pstream::parRun()` FatalError is gone.
2. **Processor patches become interfaces.** For each processor patch of phiE
   and phiI:
   - an `lduPrimitiveProcessorInterface` is built with the patch face cells
     shifted by the region offset (0 for phiE, N_E for phiI);
   - the existing `processorFvPatchField` of `phi` is passed to the solver as
     the interface field. In v2312 it reads its face cells from
     `lduAddr.patchAddr(patchId)` of the matrix it is used with, so it works on
     the combined mesh unchanged and does the MPI exchange itself;
   - `internalCoeffs` go to the diagonal, and `boundaryCoeffs` become the
     interface coefficients instead of a source term. This is exactly what
     `fvMatrix<scalar>::solveSegregated` does for coupled patches.
   Non-processor coupled patches (cyclic, AMI) still stop with a FatalError.
3. **The mesh owns the interfaces.** They are passed to the `lduPrimitiveMesh`
   constructor that takes a `PtrList<const lduInterface>` and transfers it.
   The first attempt used `addInterfaces()`, which also takes ownership, and
   that crashed with a double free when the solve returned.
4. **Log line reduced.** `pairs`, `max|eta_new - eta*|`, `I+` and `I-` are now
   summed (or maxed) over all processors.

Nothing else was needed:

- **Cell pairs are local.** `decomposeParID` puts the fluid, phiE and phiI
  copies of a cell on the same processor, so every pair link stays inside one
  processor's matrix.
- **ButlerVolmer** already reduces its total current.
- **The galvanostatic controller** (`electric::correct()`) already reduces
  `iDot`.

Message ordering: phiE and phiI may both have a processor patch towards the
same neighbour, with the same MPI tag. Both sides post their sends and
receives in the same order (phiE patches first, then phiI), and MPI does not
let messages overtake each other, so they pair up correctly.

## Running a case in parallel

Case dictionaries are the same as for the serial coupled run
(`coupled { active yes; electronRegion phiE; omega 1; }` in
`constant/phiI/regionProperties`, `phiCoupled` in `system/phiI/fvSolution`).
Keep `PCG` + `DIC` for `phiCoupled`. GAMG has not been tried on the combined
matrix.

The test case is `run/H2/HTPEMFC/implicitMethod/coupled_parallel`, a copy of
`coupled` with `parallel.sh` taken from the `phiE` case:

```bash
make mesh                                  # if the mesh is not there yet
openfoam2312 -c 'export NPROCS=4; bash ./parallel.sh' > log.parallel 2>&1
openfoam2312 -c 'mpirun -np 4 openFuelCell -parallel -fileHandler collated' > log.run 2>&1
```

## Verification

HT-PEMFC test case, 400 iterations, galvanostatic at 8000 A/m², 4 ranks
versus the serial `coupled` run (2026-10-21):

| Quantity (last iteration) | Serial | 4 ranks |
| --- | --- | --- |
| Voltage (V) | 0.4635942 | 0.4635943 |
| ibar (A/m²) | 8000.011 | 8000.011 |
| Pairs | 8000 | 8000 |
| max\|η_new − η*\| (V) | 6.480e-06 | 6.480e-06 |
| I+ / I− (A) | 0.6399991 / −0.6399991 | 0.6399991 / −0.6399991 |
| Relative imbalance | 6.2e-12 | 1.8e-12 |
| phiCoupled PCG iterations | 56 | 58 |
| Wall time (s) | 347 | 111 |

- The first coupled solve is the same in both runs: same initial residual
  (0.06062933), same currents to 6 digits.
- The voltage histories differ by at most 3.5e-4 V, early in the
  run (iteration 10: 0.48591 vs 0.48606 V). The gap shrinks to 1e-7 V by
  iteration 200. This is expected: DIC preconditioning, and the
  flow and species solvers, give slightly different iterates on 4 ranks. The
  converged answer is the same.
- No `FOAM FATAL` in `log.parallel` or `log.run`.

Not tested: restart from a written time in parallel, other rank counts, and
GAMG for `phiCoupled`.
