#!/bin/bash

set -u

# Shared with the Condor jobs; aborts if the RooUnfold build is missing.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
## Use personal setup, this one does not work currently.(I dont have RooUnfold_build_test/build/libRooUnfold.so). Might need to build RooUnfold propoerly inside workflow.
#source "${WORK}/setup_roounfold_env.sh"


# Data runs over negTagFix, so use negTagFix here too; negTag is the older production.
#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTagFix_chunks/
#INPUT_TAG=negTagFix
#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_negTag_chunks/
#INPUT_TAG=negTag
#INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/QCD/Pythia8_UParTV2_chunks/
#INPUT_TAG=UParTV2
INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/bJet/Pythia8_UParTV2_chunks/
INPUT_TAG=UParTV2
echo "Input sample ${INPUT_DIR}"


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


OUT_BASE=$mydata/bJetAggRun3/PPRef2024/bJet/agg_ntuple_chunks/
#OUT_BASE=$mydata/bJetAggRun3/PPRef2024/QCD/agg_ntuple_chunks/ ## made for tests 
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

## compilation of create_files_for_template_fit.cpp
#root -l -b -q "compile_create_files_roounfold.C(\"${ROOUNFOLD_INC}\",\"${ROOUNFOLD_BUILD}\")" > "${COMPILE_LOG}" 2>&1
##Instead:  Use direct compilation with my local setup 
root -l -b <<EOF
.L create_files_for_template_fit.cpp++
.q
EOF
# Check compilation status 
compile_status=$?
if [ ${compile_status} -ne 0 ]; then
  echo "compile failed; see ${COMPILE_LOG}"
  tail -n 80 "${COMPILE_LOG}"
  exit ${compile_status}
fi

echo "compile finished; launching chunk jobs"



#for i in $(seq 0 9); do 
## Block 9 result is missing in Matt. path. 
for i in $(seq 0 8); do  

  block=$(printf "000%d" "${i}")
  
  input="${INPUT_DIR}/merged_block_${block}_Pythia8_${INPUT_TAG}.root"
  outdir="${OUT_BASE}/block_${block}"
  ## does outdir exist? 
  # The macro picks its own filenames, so stage them and rename on the way out.
  stagedir="${outdir}/.stage${OUT_SUFFIX}"
  rm -rf "${stagedir}"
  mkdir -p "${stagedir}"

  (
  
  # Remeber: create_files_for_template_fit(Int_t RunN = 3, Int_t dataType = 2, Float_t pT_low = 80, Float_t etaCut = 2, Int_t n = 1,bool btag = true, bool isMC = true, Double_t btagWP = 0.712, bool makeTemplates = true, bool createRmatrix = true, bool makeAggNtuple = true, Long64_t ev_first = 0, Long64_t ev_last = -1, const char* inputFileOverride = "", const char* outputFolderOverride = "")
  
  #nice -n 10 root -l -b -q -e "gSystem->AddIncludePath(\"-I${ROOUNFOLD_INC} -I${ROOUNFOLD_BUILD}\"); gSystem->Load(\"${ROOUNFOLD_BUILD}/libRooUnfold.so\"); gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\")"
  # For QCD, and without default RooUnfold setup: dataType = 2
  
  #nice -n 10 root -l -b -q -e "gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,2,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\")"
  # For bJets: dataType = 1 
  nice -n 10 root -l -b -q -e "gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,1,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\")"

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
