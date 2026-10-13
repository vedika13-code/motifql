#!/bin/bash
# MotifQL final demo (Review 3). Run from the motifql folder:  bash demo/run_demo.sh
# Press Enter to move to the next step.
cd "$(dirname "$0")"
B=../build/motifql
[ -x $B ] || (cd .. && make -s)
step() { echo; echo "================ $1 ================"; read -r -p "(Enter) " _; }
search_only() { sed -n '/running searches/,$p'; }

step "1. Phases 1-2 recap: exact vs 1-mismatch search"
cat demo.motif; echo; $B demo.motif | search_only

step "2. Compiler phases, one at a time (minimised DFA shown)"
$B --dfa --no-search demo.motif | sed -n '/motif site_mm:/,$p'

step "3. Error handling"
$B bad.motif 2>&1 | grep -A2 'error:'

step "4. NEW: sensitivity - how many drifted real sites each motif still finds"
$B --sensitivity --no-search ../tests/promoter.motif 2>&1 | sed -n '/sensitivity analysis/,$p'

step "5. NEW: per-part error budgets inside a composite motif"
sed -n '/^motif/,$p' perpart.motif; echo; $B perpart.motif | search_only

step "6. NEW: the compiler chooses the budgets (budget synthesis)"
cat synth.motif; echo; $B synth.motif | sed -n '/budget synthesis/,$p'

step "7. NEW: an impossible target is rejected with a reason"
$B impossible.motif 2>&1 | head -3

step "8. NEW: errors allowed only at chosen positions"
sed -n '/^motif/,$p' positions.motif; echo; $B positions.motif | search_only

step "9. NEW: learning the budget and positions from known sites"
echo "known sites:"; grep -v '>' learn_sites.fa; echo
$B --sites learn_sites.fa --background learn_bg.fa learn.motif | sed -n '/budget synthesis/,/note:/p'

step "10. Tests and independent checker"
(cd .. && bash tests/run_tests.sh | tail -1 && python3 tools/oracle_check.py 300 7 | tail -3)
echo; echo "Demo finished."
