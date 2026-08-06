# eec_2b_analysis

CMS heavy-ion / pp-reference analysis: EEC (energy-energy correlator) between the two
b-jets in double-b events. See `README.md` for the physics-level entry points.

## Working environment: sshfs mount

Claude cannot be installed on the LLR machine, so it runs on the Mac against an sshfs
mount of that machine's home area. **The same files are visible from both sides:**

| Side | Path |
|------|------|
| Local (Mac, where Claude runs) | `/Users/zoezaidan/llruicms01/analysis_lise/eec_2b_analysis` |
| Remote (LLR, where jobs run)   | `/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis` |

The mapping is a plain prefix swap: `/Users/zoezaidan/llruicms01/...` == `/home/llr/cms/zaidan/...`.

**Do everything locally that can be done locally.** Reading, editing, writing, grepping,
and inspecting job logs all work over the mount — no SSH needed. Reach for SSH only for
commands that must execute on LLR (`condor_submit`, `condor_q`, `root`, `hadd`, anything
touching `/data_CMS` or `/cvmfs`). sshfs round-trips are slow and every SSH call pays the
ProxyJump latency, so batching remote commands into one `ssh` invocation is worth it.

Note that paths written *inside* job scripts and submit files must be the **remote**
`/home/llr/cms/zaidan/...` form, even though you are editing the file locally.

### SSH access

```bash
ssh llruicms01.in2p3.fr        # user zaidan, ProxyJump through llrgate01
```

Key-based, no password prompt; `-o BatchMode=yes` is safe. The remote `.bashrc` prints a
harmless `cmsenv: command not found` on every non-CMSSW login — filter it out of captured
output (`| grep -v cmsenv`) so it does not get mistaken for a real error.

## Runtime environment

Jobs and interactive ROOT both use the LCG view, **not** CMSSW. A CMSSW runtime that is
already loaded will collide with it, which is why the setup unsets the include and library
paths first:

```bash
unset ROOT_INCLUDE_PATH CPLUS_INCLUDE_PATH CPATH LD_LIBRARY_PATH PYTHONPATH
set +u   # LCG setup scripts are not nounset-safe
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/setup.sh
set -u
```

This gives ROOT 6.32.06. Anything that touches RooUnfold objects should instead source
`setup_roounfold_env.sh`, which does the above *and* exports `ROOUNFOLD_INC`,
`ROOUNFOLD_BUILD`, and the matching `ROOT_INCLUDE_PATH`/`LD_LIBRARY_PATH` for the private
build in `RooUnfold_build_test/`. Use it for `hadd` too — merging files that contain
RooUnfold objects fails without it.

Data lives on `/data_CMS/cms/zaidan/` (`$mydata` on the remote); shared inputs come from
`/data_CMS/cms/mnguyen/`.

## Running jobs

### HTCondor (batch, the T3 farm)

`condor_submit` and `condor_q` exist only on the remote and talk to schedd
`llrt3condor.in2p3.fr`. Submit files need the LLR T3 boilerplate:

```
Universe       = vanilla
Executable     = condor_hardprobes_data_scripts/jobs_$(Process).sh
input   = /dev/null
output  = logfiles/jobs_$(Cluster)_$(Process).sh.o
error   = logfiles/jobs_$(Cluster)_$(Process).sh.e
Log     = condor.log
getenv         = true
T3Queue        = long          # or short (the default)
WNTag          = el9

include : /opt/exp_soft/cms/t3/t3queue |

Queue 50
```

The `include : ... |` line is required — it is what injects the T3 queue configuration.

Typical cycle, generating scripts locally and submitting over SSH:

```bash
# 1. compile the shared object ONCE in the LCG env before submitting (jobs will not
#    compile it themselves; they check for it and exit 3 if missing)
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis && \
  source setup_roounfold_env.sh && \
  root -l -b -q "compile_create_files_roounfold.C(\"$ROOUNFOLD_INC\",\"$ROOUNFOLD_BUILD\")"'

# 2. generate the per-job scripts (writes condor_hardprobes_data_scripts/jobs_N.sh)
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis && \
  ./make_hardprobes_condor_scripts.sh'

# 3. submit
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis && \
  condor_submit condor.submit'

# 4. monitor
ssh llruicms01.in2p3.fr 'condor_q <cluster> -af:h JobStatus RemoteHost'
```

JobStatus: `1` idle, `2` running, `5` held. A job that disappears from `condor_q` has
finished. Then **read the logs locally** through the mount (`logfiles/`) rather than
`cat`-ing them over SSH.

Each generated job script must `cd` to the remote work dir, source
`setup_roounfold_env.sh`, `mkdir -p` its output directory, and guard against a missing
input file — see `make_hardprobes_condor_scripts.sh` for the pattern.

### Local parallel (no batch system)

`run_agg_ntuple_chunks.sh` runs the same workload as background `nice`d ROOT processes on
the interactive machine instead of Condor: it compiles once, then launches one ROOT job per
block and `wait`s. Use it for ~10-block runs; use Condor when the job count is larger.

## Conventions

- ACLiC products (`*.so`, `*.d`, `*.pcm`), ROOT files, plots, and logs are gitignored.
  Do not commit them.
- Compile the shared object once and let jobs `gSystem->Load()` it. Parallel ACLiC
  compilation in the same directory produces half-written artifacts; the run scripts
  delete stale ones before recompiling for this reason.
- Emacs backup and lock files (`*~`, `#*#`) litter the tree and are gitignored — ignore them.
- `._*` files are macOS AppleDouble metadata created by the mount, not real source files.
