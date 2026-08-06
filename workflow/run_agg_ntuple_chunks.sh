#!/bin/bash

set -u

# Shared with the Condor jobs; aborts if the RooUnfold build is missing.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
source "${WORK}/setup_roounfold_env.sh"


# Data runs over negTagFix, so use negTagFix here too; negTag is the older production.
#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTagFix_chunks/
#INPUT_TAG=negTagFix
#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTag_chunks/
#INPUT_TAG=negTag
INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_UParTV2_chunks/
INPUT_TAG=UParTV2


# Set OUT_TAG= (empty) for the original untagged filenames.
OUT_TAG=upartv2
OUT_SUFFIX="${OUT_TAG:+_${OUT_TAG}}"

# Must match BTAG_WP in make_hardprobes_condor_scripts.sh.
case "${INPUT_TAG}" in
  UParTV2)   BTAG_WP=0.712 ;;
  negTagFix) BTAG_WP=0.868 ;;
  negTag)    BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac
echo "sample ${INPUT_TAG}, b-tag WP ${BTAG_WP}"

OUT_BASE=$mydata/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/
LOG_DIR=${OUT_BASE}/logs

# Must match the build dir set in rootlogon.C.
ACLIC_BUILD_DIR=$mydata/bJetAggRun3/PPRef2024/build

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Drop half-written ACLiC products from a previous parallel launch.
rm -f "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.d \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.so \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_dict.* \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_map.*

root -l -b -q "compile_create_files_roounfold.C(\"${ROOUNFOLD_INC}\",\"${ROOUNFOLD_BUILD}\")" > "${COMPILE_LOG}" 2>&1
compile_status=$?
if [ ${compile_status} -ne 0 ]; then
  echo "compile failed; see ${COMPILE_LOG}"
  tail -n 80 "${COMPILE_LOG}"
  exit ${compile_status}
fi

echo "compile finished; launching chunk jobs"



for i in $(seq 0 9); do  
  block=$(printf "000%d" "${i}")
  
  input="${INPUT_DIR}/merged_block_${block}_Pythia8_${INPUT_TAG}.root"
  outdir="${OUT_BASE}/block_${block}"
  # The macro picks its own filenames, so stage them and rename on the way out.
  stagedir="${outdir}/.stage${OUT_SUFFIX}"
  rm -rf "${stagedir}"
  mkdir -p "${stagedir}"

  (
    nice -n 10 root -l -b -q -e "gSystem->AddIncludePath(\"-I${ROOUNFOLD_INC} -I${ROOUNFOLD_BUILD}\"); gSystem->Load(\"${ROOUNFOLD_BUILD}/libRooUnfold.so\"); gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\")"

    for f in "${stagedir}"/*.root; do
      [ -e "${f}" ] || continue
      base=$(basename "${f}")
      mv "${f}" "${outdir}/${base%.root}${OUT_SUFFIX}.root"
    done
    rmdir "${stagedir}" 2>/dev/null
  ) > "${LOG_DIR}/block_${block}${OUT_SUFFIX}.log" 2>&1 &


  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
