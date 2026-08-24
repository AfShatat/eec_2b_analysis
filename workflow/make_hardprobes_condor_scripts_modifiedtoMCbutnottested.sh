#!/bin/bash

set -u

# WORK=/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow
WORK=/home/llr/cms/shatat/CMSAnalysis/eec_2b_analysis/workflow ## bJet
# OUT_BASE=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks
OUT_BASE=/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/bJet ## bJet


# ACLIC_BUILD_DIR must match the build dir set in rootlogon.C.
# ACLIC_BUILD_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/build
ACLIC_BUILD_DIR=/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/build
# SCRIPT_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/condor/hardprobes_data_scripts
SCRIPT_DIR=/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/condor/bJet_condor_scripts
# CONDOR_LOG_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/condor/logfiles
CONDOR_LOG_DIR=/data_CMS/cms/shatat/bJetAggRun3/PPRef2024/condor/logfiles

# INPUT_TAG names the chunk subdir (<TAG>_chunks) and the per-block filenames.
#INPUT_TAG=negTagFix
INPUT_TAG=UParTV2

# Suffix on every output .root; leave empty for untagged names.
OUT_TAG=upartv2
OUT_SUFFIX="${OUT_TAG:+_${OUT_TAG}}"

# Must match BTAG_WP in run_agg_ntuple_chunks.sh.
#  UParTV2)    BTAG_WP=0.872 ;; // Like Zoe templates 
#   UParTV2)    BTAG_WP=0.712 ;; // git version --- latest and correct. 
case "${INPUT_TAG}" in
  UParTV2)    BTAG_WP=0.712 ;; 
  negTagFix)  BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac
echo "sample ${INPUT_TAG}, b-tag WP ${BTAG_WP}"

mkdir -p "${SCRIPT_DIR}" "${CONDOR_LOG_DIR}"
rm -f "${SCRIPT_DIR}"/job_*.sh "${SCRIPT_DIR}"/jobs_*.sh

## Change loop for MC 
job_index=0
#for primary_dataset in $(seq 0 4); do
#  input_dir=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/${primary_dataset}/${INPUT_TAG}_chunks
for i in $(seq 0 9); do
  block=$(printf "000%d" "${i}")
  script="${SCRIPT_DIR}/jobs_${job_index}.sh"
  input="${input_dir}/merged_block_${block}_${INPUT_TAG}.root"
  outdir="${OUT_BASE}/HardProbes${primary_dataset}/block_${block}"
  logdir="${OUT_BASE}/HardProbes${primary_dataset}/logs"
  # The macro picks its own filenames, so stage them and rename on the way out.
  stagedir="${outdir}/.stage${OUT_SUFFIX}"

  cat > "${script}" <<EOF

#!/bin/bash

set -euo pipefail

WORK=${WORK}
JOB_INDEX=${job_index}

# PRIMARY_DATASET=${primary_dataset}

BLOCK=${block}
INPUT=${input}
OUTDIR=${outdir}
LOGDIR=${logdir}
STAGEDIR=${stagedir}
OUT_SUFFIX=${OUT_SUFFIX}
ACLIC_BUILD_DIR=${ACLIC_BUILD_DIR}

cd "\${WORK}"
source "\${WORK}/setup_roounfold_env.sh"
mkdir -p "\${OUTDIR}" "\${LOGDIR}"
rm -rf "\${STAGEDIR}"
mkdir -p "\${STAGEDIR}"

if [ ! -f "\${INPUT}" ]; then
  echo "ERROR: missing input \${INPUT}"
  exit 2
fi

if [ ! -f "\${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so" ]; then
  echo "ERROR: missing \${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so -- run ./run_agg_ntuple_chunks.sh style compile (LCG env) once before submitting"
  exit 3
fi

root -l -b -q -e "gSystem->AddIncludePath(\"-I\${ROOUNFOLD_INC} -I\${ROOUNFOLD_BUILD}\"); if (gSystem->Load(\"\${ROOUNFOLD_BUILD}/libRooUnfold.so\") < 0) gSystem->Exit(3); if (gSystem->Load(\"\${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\") < 0) gSystem->Exit(3); create_files_for_template_fit(3,2,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"\${INPUT}\",\"\${STAGEDIR}\")"

for f in "\${STAGEDIR}"/*.root; do
  [ -e "\${f}" ] || continue
  base=\$(basename "\${f}")
  mv "\${f}" "\${OUTDIR}/\${base%.root}\${OUT_SUFFIX}.root"
done
rmdir "\${STAGEDIR}" 2>/dev/null || true
EOF

  chmod +x "${script}"
  job_index=$((job_index + 1))
done

echo "created 50 scripts in ${SCRIPT_DIR}"


## Check now the compialtion?? and run_agg.sh