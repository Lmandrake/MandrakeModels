"""Run a .NET Framework exe on RimWorld's own Mono (the game's embedded mono-2.0-bdwgc.dll and its Managed class
library), the way Unity embeds it. Windows Python only:

    python.exe monohost.py LiTiming.exe <args...>

Only reads the game folder; nothing in it is written. Exit code is the program's.
"""
import ctypes
import os
import sys

RIM = os.environ.get("RIMWORLD_DIR", r"C:\Program Files (x86)\Steam\steamapps\common\RimWorld")
DLL = os.path.join(RIM, r"MonoBleedingEdge\EmbedRuntime\mono-2.0-bdwgc.dll")
MANAGED = os.path.join(RIM, r"RimWorldWin64_Data\Managed")
ETC = os.path.join(RIM, r"MonoBleedingEdge\etc")


def main():
    exe = os.path.abspath(sys.argv[1])
    m = ctypes.CDLL(DLL)
    vp, cp = ctypes.c_void_p, ctypes.c_char_p
    m.mono_set_dirs.argtypes = [cp, cp]
    m.mono_set_assemblies_path.argtypes = [cp]
    m.mono_jit_init_version.argtypes = [cp, cp]
    m.mono_jit_init_version.restype = vp
    m.mono_domain_assembly_open.argtypes = [vp, cp]
    m.mono_domain_assembly_open.restype = vp
    m.mono_jit_exec.argtypes = [vp, vp, ctypes.c_int, ctypes.POINTER(cp)]
    m.mono_jit_exec.restype = ctypes.c_int
    m.mono_set_dirs(MANAGED.encode(), ETC.encode())
    m.mono_set_assemblies_path(MANAGED.encode())
    m.mono_config_parse(None)
    dom = m.mono_jit_init_version(b"LiTiming", b"v4.0.30319")
    if not dom:
        sys.exit("mono_jit_init_version failed")
    asm = m.mono_domain_assembly_open(dom, exe.encode())
    if not asm:
        sys.exit("could not open " + exe)
    args = [exe] + sys.argv[2:]
    argv = (cp * len(args))(*[a.encode() for a in args])
    sys.stdout.flush()
    rc = m.mono_jit_exec(dom, asm, len(args), argv)
    os._exit(rc)


if __name__ == "__main__":
    main()
