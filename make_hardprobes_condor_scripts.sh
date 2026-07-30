#!/bin/bash

set -u

WORK=/home/llr/cms/mnguyen/eec_2b_analysis
SCRIPT_DIR=${WORK}/condor_hardprobes_data_scripts
OUT_BASE=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks

mkdir -p "${SCRIPT_DIR}"
rm -f "${SCRIPT_DIR}"/job_*.sh "${SCRIPT_DIR}"/jobs_*.sh

job_index=0
for primary_dataset in $(seq 0 4); do
  input_dir=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/${primary_dataset}/negTagFix_chunks

  for i in $(seq 0 9); do
    block=$(printf "000%d" "${i}")
    script="${SCRIPT_DIR}/jobs_${job_index}.sh"
    input="${input_dir}/merged_block_${block}_negTagFix.root"
    outdir="${OUT_BASE}/HardProbes${primary_dataset}/block_${block}"
    logdir="${OUT_BASE}/HardProbes${primary_dataset}/logs"

    cat > "${script}" <<EOF
#!/bin/bash

set -euo pipefail

WORK=${WORK}
JOB_INDEX=${job_index}
PRIMARY_DATASET=${primary_dataset}
BLOCK=${block}
INPUT=${input}
OUTDIR=${outdir}
LOGDIR=${logdir}
LOCK=\${WORK}/create_files_for_template_fit.compile.lock

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /grid_mnt/vol_home/llr/cms/mnguyen/prod/forest/bJetAgg2024PPRef/CMSSW_14_1_9_patch2/src
eval \$(scram runtime -sh)

cd "\${WORK}"
mkdir -p "\${OUTDIR}" "\${LOGDIR}"

if [ ! -f "\${INPUT}" ]; then
  echo "ERROR: missing input \${INPUT}"
  exit 2
fi

(
  flock 9
  root -l -b -q -e '.L create_files_for_template_fit.cpp+' > "\${LOGDIR}/compile_jobs_${job_index}.log" 2>&1
) 9>"\${LOCK}"

root -l -b -q -e "gSystem->Load(\"\${WORK}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,0,80,200,2,1,true,false,0.868,true,false,true,0,-1,\"\${INPUT}\",\"\${OUTDIR}\")"
EOF

    chmod +x "${script}"
    job_index=$((job_index + 1))
  done
done

echo "created 50 scripts in ${SCRIPT_DIR}"
