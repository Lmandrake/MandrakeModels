# The C# port (stage C)

`csharp/` holds a standalone console port of `local_inertial` and `local_inertial_frcap`
(`src/mandrakemodels/solvers/local_inertial.py`, with `base.py`'s transfer and forcing). Its purpose was to
measure the solver's cost on the runtime RimWorld uses: Unity's embedded Mono. It is not mod code, and nothing
in it writes to the game folder.

## What is in it

| file | role |
|---|---|
| `csharp/LiTiming/LocalInertialDL.cs` | the solver, on integer water |
| `csharp/LiTiming/SceneFile.cs` | reads the `.lis` scene files written by `python -m mandrakemodels.csharp_export` |
| `csharp/LiTiming/Spreader.cs` | spreads one recompute over its 12 game ticks (weighted block passes), bitwise equal to running it whole |
| `csharp/LiTiming/Program.cs` | the console: subcommands below |
| `csharp/LiTiming/LiTiming.csproj` | builds `net8.0` and `net472` |
| `csharp/MonoHost/monohost.py` | runs the `net472` build on RimWorld's own Mono, embedded from Windows Python by ctypes the way Unity embeds it |
| `csharp/run_stagec.sh` | the whole run: export, build, every subcommand on three runtimes |
| `csharp/summarize_stagec.py` | summarises the run's output |

### The solver on integer water

- Water is an `int` per square, in decilitres (dL).
- Each face keeps a `double` discharge (the momentum memory) and a sub-dL rounding residue, carried to the next
  step.
- The donor limiter works in integers, and every face moves whole dL once with opposite signs, so the books
  are exact to the dL.
- Wet-only stepping works in 16 x 16 blocks, dilated by one block, and gives bitwise the same answer as the
  full grid.
- An optional sleep mode lets settled water stop being stepped.

## Scene files

`python -m mandrakemodels.csharp_export --out DIR [--sizes 128,250] [--scenes ...] [--solver NAME] [--no-ref]`
writes `DIR/<scene>_<n>.lis`. This is a little-endian binary file, magic `LIS2`, holding:

- the grid, Manning n, theta, alpha, the reference solver's face Froude cap (0 = none) and the duration
- the bed, depth and walls
- the forcing
- the dig events
- unless `--no-ref`: the Python reference's depth at 4 evenly spaced checkpoints, with step counts and volume
  error

`--solver` is one of `local_inertial`, `local_inertial_frcap` or `local_inertial_frcap07`. The only callable
forcing it can encode is the swell sine. Any other callable forcing, and rain or groundwater in forms it
cannot encode, is refused with an error. Besides the standard scenes it knows two timing loads,
`census_river` and `stress_half` (see `scenes.md`).

## Subcommands

`LiTiming <cmd> ...` (from the doc comment in `Program.cs`):

| command | does |
|---|---|
| `fidelity DIR OUT.jsonl` | every scene against the Python checkpoints |
| `timing DIR OUT.csv [--modes full,wet,sleep] [--only a,b] [--label L]` | cost per game tick: a recompute every 12 ticks = 0.2 water-s, in CFL substeps |
| `stability DIR OUT.jsonl` | long runs |
| `fuzz FILE.lis OUT.jsonl [--seqs 1000]` | random digs, dams, walls and demolitions; the books must stay exact |
| `spread DIR OUT.csv` | recompute spread over 12 ticks, checked bitwise against whole |
| `settled DIR OUT.csv` | cost once water sleeps |
| `wake FILE.lis OUT.jsonl` | a dig into sleeping water must flow as in water that never slept |
| `evensteps FILE.lis` | every driver takes equal steps (the rule at `LocalInertialDL.EvenDt`) |
| `profile`, `drift` | diagnostics |

`LITIMING_MUTATE=manning|theta|walldelete|frcap1|nocap|nowake|nodeadband` plants a known fault, to prove
that the checks go red. `nocbrt` is a timing probe.

## Running it

`run_stagec.sh` was written for one setup: WSL on Windows, building with a Windows `dotnet.exe`, and running
Mono from a Windows Python against a local RimWorld install. Its defaults name paths on the author's machine
(`STAGEC=/mnt/d/...`, `DOTNET=/mnt/c/Users/...`). To use it elsewhere, set:

- `STAGEC`: the work and output directory
- `DOTNET`: the dotnet executable
- `RIMWORLD_DIR`: for `monohost.py`; its default is the usual Steam install path
- `REF`: the Python reference solver (default `local_inertial_frcap`)
- `STEPS`: a subset of `evensteps,fidelity,timing,spread,settled,wake,stability,fuzz`
- `SIZES`: the timing sizes

It also calls `wslpath` and `python.exe`, so outside WSL it needs editing. The `.NET 8` build alone should
run anywhere `dotnet build` works, as `dotnet LiTiming.dll <cmd> ...`. This was not tested for this document:
the machine it was written on has no `dotnet`.

## Results

`results/stagec.md`: fidelity (40 of 40 rows within tolerance on Unity Mono and .NET 8), cost per tick by map
size and runtime, the bars it was held to, and the findings, including why every driver must take equal
substeps. Those figures were measured on the author's Windows machine. The raw output paths it cites are on
that machine and are not in this repository.
