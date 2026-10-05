#!/usr/bin/env bash
# =============================================================================
# run_leave_one_out.sh
#
# Runs the four leave-one-out detections for the AKOYA circularity audit, then
# restores LEAVE_ONE_OUT = None.
#
# WHY A DRIVER AND NOT FOUR MANUAL EDITS
#     The one real hazard in this procedure is leaving LEAVE_ONE_OUT set and
#     rerunning script 04 later. That would write a sensitivity run into
#     structures_rev5 and every downstream script would silently inherit it.
#     The restore here is on an EXIT trap, so it runs on success, on failure,
#     and on Ctrl-C.
#
# WHAT IT DOES NOT DO
#     It does not touch structures_rev5. Verified: in rev6 the OUT_DIR redirect
#     is at line 394 and FIG_DIR / TAB_DIR / CELL_DIR are built from it at lines
#     408 to 410, with no file written before the redirect. So a leave-one-out
#     run cannot reach the primary directory.
#
# USAGE
#     cd /master/jlehle/WORKING/AKOYA/scripts
#     conda activate sc_pre
#     bash run_leave_one_out.sh
#     python AKOYA_04g_Leave_One_Out_Audit.py
#
#     Each detection takes about as long as a normal script 04 run, and there
#     are four, so expect roughly four times that. Nothing runs in parallel on
#     purpose: four concurrent runs on one node is how a job gets killed for
#     memory halfway through and leaves a half-written directory.
# =============================================================================
set -euo pipefail

SCRIPT="AKOYA_04_Structures_rev6.py"
BACKUP="${SCRIPT}.loo_backup"

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: $SCRIPT not found. Run this from the scripts directory."
    exit 1
fi

# Refuse to start from a file that is already set, because then the backup
# would restore the wrong value and the primary run would stay redirected.
CURRENT=$(grep -m1 '^LEAVE_ONE_OUT' "$SCRIPT" || true)
if [[ "$CURRENT" != "LEAVE_ONE_OUT = None" ]]; then
    echo "ERROR: expected 'LEAVE_ONE_OUT = None' in $SCRIPT, found:"
    echo "       $CURRENT"
    echo "       Set it back to None before running this driver, so the"
    echo "       restore at the end puts back the right value."
    exit 1
fi

cp "$SCRIPT" "$BACKUP"
restore() {
    if [[ -f "$BACKUP" ]]; then
        cp "$BACKUP" "$SCRIPT"
        rm -f "$BACKUP"
        echo ""
        echo "=== restored: $(grep -m1 '^LEAVE_ONE_OUT' "$SCRIPT")"
    fi
}
trap restore EXIT

PHENOS=(
    "CD68+IDO1+ Macrophages"
    "CD68+IDO1- Macrophages"
    "CD163+ Macrophages"
    "Neutrophils"
)

for P in "${PHENOS[@]}"; do
    echo ""
    echo "============================================================"
    echo "LEAVE-ONE-OUT: $P"
    echo "============================================================"
    cp "$BACKUP" "$SCRIPT"
    # python rather than sed: the phenotype names contain + and spaces, and the
    # assert makes a silent no-op substitution impossible.
    python - "$SCRIPT" "$P" <<'PY'
import sys, re
path, pheno = sys.argv[1], sys.argv[2]
src = open(path).read()
new, n = re.subn(r'(?m)^LEAVE_ONE_OUT\s*=.*$',
                 'LEAVE_ONE_OUT = ' + repr(pheno), src)
if n != 1:
    raise SystemExit(f"ERROR: expected exactly one LEAVE_ONE_OUT assignment, "
                     f"found {n}. Not writing.")
open(path, 'w').write(new)
print(f"    set LEAVE_ONE_OUT = {pheno!r}")
PY
    python "$SCRIPT"
done

echo ""
echo "============================================================"
echo "All four detections done. Directories written:"
ls -d ../structures_rev5_no_* 2>/dev/null || \
    ls -d /master/jlehle/WORKING/AKOYA/structures_rev5_no_* 2>/dev/null || \
    echo "    none found, which means the runs did not write where expected"
echo ""
echo "Next:  python AKOYA_04g_Leave_One_Out_Audit.py"
echo "============================================================"
