#!/bin/bash

set -u

WORK=/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis
SCRIPT_DIR=${WORK}/condor_hardprobes_data_scripts
OUT_BASE=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks

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

cd "\${WORK}"
source "\${WORK}/setup_roounfold_env.sh"
mkdir -p "\${OUTDIR}" "\${LOGDIR}"

if [ ! -f "\${INPUT}" ]; then
  echo "ERROR: missing input \${INPUT}"
  exit 2
fi

if [ ! -f "\${WORK}/create_files_for_template_fit_cpp.so" ]; then
  echo "ERROR: missing \${WORK}/create_files_for_template_fit_cpp.so -- run ./run_agg_ntuple_chunks.sh style compile (LCG env) once before submitting"
  exit 3
fi

root -l -b -q -e "gSystem->AddIncludePath(\"-I\${ROOUNFOLD_INC} -I\${ROOUNFOLD_BUILD}\"); if (gSystem->Load(\"\${ROOUNFOLD_BUILD}/libRooUnfold.so\") < 0) gSystem->Exit(3); if (gSystem->Load(\"\${WORK}/create_files_for_template_fit_cpp.so\") < 0) gSystem->Exit(3); create_files_for_template_fit(3,0,80,120,2,1,true,false,0.868,true,false,true,0,-1,\"\${INPUT}\",\"\${OUTDIR}\")"
EOF

    chmod +x "${script}"
    job_index=$((job_index + 1))
  done
done

echo "created 50 scripts in ${SCRIPT_DIR}"
