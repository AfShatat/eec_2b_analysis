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


INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTag_chunks/
OUT_BASE=$mydata/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/
LOG_DIR=${OUT_BASE}/logs

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Avoid stale or half-written ACLiC products from a previous parallel launch.
rm -f create_files_for_template_fit_cpp.d \
      create_files_for_template_fit_cpp.so \
      create_files_for_template_fit_cpp_ACLiC_dict.* \
      create_files_for_template_fit_cpp_ACLiC_map.*

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
  
  input="${INPUT_DIR}/merged_block_${block}_Pythia8_negTag.root"
  outdir="${OUT_BASE}/block_${block}"
  mkdir -p "${outdir}"

  nice -n 10 root -l -b -q -e "gSystem->AddIncludePath(\"-I${ROOUNFOLD_INC} -I${ROOUNFOLD_BUILD}\"); gSystem->Load(\"${ROOUNFOLD_BUILD}/libRooUnfold.so\"); gSystem->Load(\"${WORK}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,120,2,1,true,true,0.868,true,true,true,0,-1,\"${input}\",\"${outdir}\")" \
   > "${LOG_DIR}/block_${block}.log" 2>&1 &
       

  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
