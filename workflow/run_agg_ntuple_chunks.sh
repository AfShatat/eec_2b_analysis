#!/bin/bash

#set -u

# RooUnfold setup. Keep this self-contained so batch shells match interactive ROOT.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
ROOUNFOLD_DIR=${WORK}/RooUnfold_build_test
ROOUNFOLD_INC=${ROOUNFOLD_DIR}/src/src
ROOUNFOLD_BUILD=${ROOUNFOLD_DIR}/build

unset ROOT_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH
unset CPATH
# Avoid mixing an already-loaded CMSSW runtime with the LCG view.
unset LD_LIBRARY_PATH
unset PYTHONPATH

# LCG setup scripts are not nounset-safe.
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/setup.sh
set -u

export ROOT_INCLUDE_PATH=${ROOUNFOLD_INC}:${ROOUNFOLD_BUILD}
export LD_LIBRARY_PATH=${ROOUNFOLD_BUILD}:${LD_LIBRARY_PATH:-}


#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTag_chunks/
#INPUT_TAG=negTag
INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_UParTV2_chunks/
INPUT_TAG=UParTV2


# To switch back to negTag: uncomment the two negTag lines above, comment the two
# UParTV2 ones, and set OUT_TAG= (empty) to get the original untagged filenames
OUT_TAG=upartv2
OUT_SUFFIX="${OUT_TAG:+_${OUT_TAG}}"

# b-tag working point follows the tagger version, so it cannot drift out of sync
# with the sample selected above: UParT v2 -> 0.872, v1 (negTag) -> 0.868.
case "${INPUT_TAG}" in
  UParTV2) BTAG_WP=0.872 ;;
  negTag)  BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac
echo "sample ${INPUT_TAG}, b-tag WP ${BTAG_WP}"

OUT_BASE=$mydata/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/
LOG_DIR=${OUT_BASE}/logs

# ACLiC build products live on /data_CMS, not next to the code. Must match the
# path set in rootlogon.C, which is what actually redirects them.
ACLIC_BUILD_DIR=$mydata/bJetAggRun3/PPRef2024/build

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Avoid stale or half-written ACLiC products from a previous parallel launch.
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
  # create_files_for_template_fit picks its own output filenames, so let it write
  # into a private staging dir and rename on the way out.
  stagedir="${outdir}/.stage${OUT_SUFFIX}"
  rm -rf "${stagedir}"
  mkdir -p "${stagedir}"

  (
    nice -n 10 root -l -b -q -e "gSystem->AddIncludePath(\"-I${ROOUNFOLD_INC} -I${ROOUNFOLD_BUILD}\"); gSystem->Load(\"${ROOUNFOLD_BUILD}/libRooUnfold.so\"); gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,120,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\")"

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
