import json
import os
from pathlib import Path
import sys
import subprocess

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import configure_belly as c


def touch_executable(path):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)


def test_config_uses_absolute_prefix_and_existing_canonical_envs(tmp_path):
    envs=tmp_path/"envs"
    for tool,(names,executables) in c.SPECS.items():
        touch_executable(envs/names[0]/"bin"/executables[0])
    touch_executable(envs/"mito-novoplasty/bin/NOVOPlasty4.3.5.pl")
    itsme=tmp_path/"itsme.v2.sh"
    touch_executable(itsme)
    conda=tmp_path/"bin/conda"
    touch_executable(conda)
    config,rows=c.build_config(envs,conda,itsme)
    assert all(r["status"]=="found" for r in rows)
    for tool,cmd in config["commands"].items():
        assert cmd[:4]==[str(conda),"run","--no-capture-output","-p"]
        assert str(envs) in cmd[4]
    assert config["commands"]["getorganelle"][4]==str(envs/"getorganelle")
    assert config["commands"]["metaspades"][4]==str(envs/"itsme")


def test_novoplasty_is_located_in_conda_package_share_directory(tmp_path):
    prefix=tmp_path/"mito-novoplasty"
    path=prefix/"share/novoplasty-4.3.5/NOVOPlasty4.3.5.pl"
    path.parent.mkdir(parents=True)
    path.write_text("# Perl script\n")
    touch_executable(prefix/"bin/perl")
    meta=prefix/"conda-meta/novoplasty-4.3.5-0.json"
    meta.parent.mkdir()
    meta.write_text(json.dumps({"name":"novoplasty","files":[str(path.relative_to(prefix))]}))
    assert c.novo_executable(prefix)==[str(prefix/"bin/perl"),str(path)]


def test_resolution_falls_back_to_existing_alternative_and_reports_missing(tmp_path):
    touch_executable(tmp_path/"getorganelle/bin/get_organelle_from_reads.py")
    selected,alternatives=c.locate(tmp_path,["mito-getorganelle","getorganelle"],["get_organelle_from_reads.py"])
    assert selected[0]==tmp_path/"getorganelle"
    assert alternatives==[str(tmp_path/"getorganelle")]
    assert c.locate(tmp_path,["nonexistent"],["mitofinder"])[0] is None


def test_reference_scan_is_bounded_and_ignores_environment_trees(tmp_path):
    ref=tmp_path/"mitos/refseq89m"
    ref.mkdir(parents=True)
    (ref/"auxinfo.json").write_text("{}")
    ignored=tmp_path/"envs/fake/refseq39"
    ignored.mkdir(parents=True)
    (ignored/"auxinfo.json").write_text("{}")
    report=c.reference_inventory([tmp_path])
    paths=[x["path"] for x in report["candidates"]]
    assert str(ref) in paths and str(ignored) not in paths
    assert c.reference_inventory([tmp_path],max_dirs=1)["truncated"]


def conda_double(path):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text('''#!/bin/sh
test "$1" = run || exit 90
shift
test "$1" = --no-capture-output && shift
test "$1" = -p || exit 91
shift
CONDA_PREFIX="$1"
export CONDA_PREFIX
PATH="$CONDA_PREFIX/bin:$PATH"
export PATH
shift
if test -n "${MITO_TEST_CONDA_TRACE:-}"; then
  printf '%s\\n' "$CONDA_PREFIX" >> "$MITO_TEST_CONDA_TRACE"
fi
exec "$@"
''')
    path.chmod(0o755)


def test_runtime_switches_each_stage_and_records_first_choice_alternatives(tmp_path):
    envs=tmp_path/"envs with spaces"
    for env,tool in [("itsme","fastp"),("getorganelle","get_organelle_from_reads.py"),
                     ("mito-getorganelle","get_organelle_from_reads.py"),
                     ("mito-mtgrasp","mtgrasp.py"),("mtgrasp","mtgrasp.py")]:
        p=envs/env/"bin"/tool
        touch_executable(p)
        p.write_text('#!/bin/sh\nprintf "%s\\n" "$CONDA_PREFIX"\n')
    conda=tmp_path/"bin/conda"
    conda_double(conda)
    commands,rows=c.runtime_commands(envs,conda,tmp_path/"missing.sh",["fastp","getorganelle","mtgrasp"])
    expected=[envs/"itsme",envs/"getorganelle",envs/"mito-mtgrasp"]
    for tool,prefix in zip(["fastp","getorganelle","mtgrasp"],expected):
        result=subprocess.run(commands[tool],capture_output=True,text=True,check=True)
        assert result.stdout.strip()==str(prefix)
    getorg=next(r for r in rows if r["tool"]=="getorganelle")
    assert str(envs/"mito-getorganelle") in getorg["note"]


def test_launcher_enters_controller_environment_without_manual_activation(tmp_path):
    root=Path(__file__).resolve().parents[1]
    envs=tmp_path/"envs"
    (envs/"itsme/bin").mkdir(parents=True)
    (envs/"itsme/bin/python").symlink_to(sys.executable)
    conda_double(tmp_path/"bin/conda")
    trace=tmp_path/"trace.txt"
    env=dict(os.environ,MITO_ENV_ROOT=str(envs),MITO_TEST_CONDA_TRACE=str(trace))
    env.pop("MITO_PYTHON",None)
    env.pop("MITO_CONDA_EXE",None)
    p=subprocess.run(["bash",str(root/"mito_looper.sh"),"--version"],env=env,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    assert p.stdout.strip()=="2.0.8"
    assert trace.read_text().strip()==str(envs/"itsme")
