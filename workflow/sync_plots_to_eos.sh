#!/bin/bash

# Copy plots from /data_CMS to the CERNBox website on EOS. Run on llruicms01.
#
#   ./sync_plots_to_eos.sh                      # upartv2 template fit -> template_fit_UParTv2
#   ./sync_plots_to_eos.sh <src_dir> <name>     # any other plot folder
#
# Needs a Kerberos ticket: kinit zzaidanc@CERN.CH

set -u

CERN_USER=zzaidanc
EOS_WWW=/eos/user/${CERN_USER:0:1}/${CERN_USER}/www

RESULTS=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results

SRC=${1:-${RESULTS}/TemplateFits_Run3_minHLT60_LinearBin_upartv2_www}
DEST_NAME=${2:-template_fit_UParTv2}
DEST=${EOS_WWW}/${DEST_NAME}

if [ ! -d "${SRC}" ]; then
  echo "no such source folder: ${SRC}"
  echo "(run this on llruicms01 -- /data_CMS is not visible from the Mac mount)"
  exit 1
fi

n_png=$(ls "${SRC}"/*.png 2>/dev/null | wc -l)
if [ "${n_png}" -eq 0 ]; then
  echo "no pngs in ${SRC}; nothing to copy"
  exit 1
fi

# A missing ticket only means a password prompt, so warn rather than fail.
if ! klist -s 2>/dev/null; then
  echo "warning: no Kerberos ticket (run 'kinit ${CERN_USER}@CERN.CH' to avoid password prompts)"
fi

echo "copying ${n_png} pngs"
echo "  from ${SRC}"
echo "  to   ${CERN_USER}@lxplus.cern.ch:${DEST}"

ssh "${CERN_USER}@lxplus.cern.ch" "mkdir -p ${DEST}" || exit 1

# --no-perms/--no-group: EOS fuse rejects rsync's ownership calls.
rsync -rltvz --no-perms --no-group \
  "${SRC}/" "${CERN_USER}@lxplus.cern.ch:${DEST}/" || exit 1

echo
echo "done: https://${CERN_USER}.web.cern.ch/${DEST_NAME}/"
echo "(directory listing is off by default -- add CERN's index.php to browse it as a gallery)"
