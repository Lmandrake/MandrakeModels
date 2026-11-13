# Stage C: local-inertial in C#, on integer decilitres

A standalone console port of `src/mandrakemodels/solvers/local_inertial.py` (with `base.py`'s transfer and forcing),
for timing on the runtime RimWorld uses. Not mod code; nothing here touches the game folder except reading its Mono.

- `LiTiming/LocalInertialDL.cs`: the solver. Water is an int per square in dL; each face keeps a double discharge
  (the momentum memory) and a sub-dL rounding residue carried to the next step; integer donor limiter; every face moves
  whole dL once with opposite signs, so the books are exact. Wet-only stepping in 16x16 blocks (dilated by one block)
  gives bitwise the same answer as the full grid.
- `LiTiming/Program.cs`: `fidelity` (against Python checkpoints), `timing` (game clock: a recompute every 12 ticks =
  0.2 water-s, CFL substeps), `stability`, `fuzz`, `profile`.
- `MonoHost/monohost.py`: runs the net472 build on **RimWorld's own Mono** (`MonoBleedingEdge\EmbedRuntime\mono-2.0-bdwgc.dll`
  with `RimWorldWin64_Data\Managed` as the class library), embedded from Windows Python by ctypes the way Unity embeds it.
- `run_stagec.sh`: the whole run (build on `D:\ClaudeWork\mandrakemodels\stagec`, since dotnet.exe cannot build from a
  WSL path).

Scenes:

    python -m mandrakemodels.csharp_export --out $STAGEC/scenes --sizes 128,250               # with Python checkpoints
    python -m mandrakemodels.csharp_export --out $STAGEC/load --sizes 250,325,400 --no-ref \
        --scenes river,lake_at_rest,pit_dig,dam_break,swell,census_river,stress_half          # timing loads

`LITIMING_MUTATE=manning|theta|walldelete` plants a fault (for proving the checkers go red); `nocbrt` is a timing probe.
Results: `results/stagec.md`.
