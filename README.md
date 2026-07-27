# Calculate PRS metrics for disease-progression trajectories

This is a shareable pipeline for evaluating polygenic scores derived from **time-to-event GWAS**
(e.g. GATE) across clinically ordered disease trajectories `T0 → T1 → T2`. 

Example — the `T2D_to_CAD` trajectory: T1 = type-2 diabetes, T2 = coronary artery disease.

This pipeline uses **three** PRS, one for each step of the sequence:

- **onset PRS** (`PRS_T0T1`) — genetic risk of reaching T1;
- **outcome PRS** (`PRS_T0T2`) — genetic risk of the final condition T2 in the general population;
- **progression PRS** (`PRS_T1T2`) — genetic risk of moving from T1 to T2, learned among people who
  already have T1.

**The question this answers** (two sides of it): among people who already have T1, does the
*progression* PRS add information about who will go on to T2 — **beyond the onset PRS** and **beyond
the general outcome PRS**? To answer it, the pipeline fits a ladder of eight models (M0–M7, below)
that add the three scores in different combinations and compares them. The two matching comparisons —
**M4-M1** (progression beyond onset) and **M5-M2** (progression beyond outcome) — are the two
**co-primary** tests; **M7−M6** (progression beyond onset *and* outcome) is the stronger fully-adjusted
confirmation.

**Two separate analyses ("layers"), each run as its own job:**

- **Layer 1 — association (GWAS-aligned).** Looks back over recorded T1→T2 times among people with
  T1. Asks: *is genetic burden associated with a shorter recorded time from T1 to T2?*
- **Layer 2 — prospective prediction.** Starts the clock at a defined index (recruitment, or the T1
  date) among people who are T2-free at that point. Asks: *does the score predict future T2?*
  
  
  Note **Layer 2** analysis requires information regarding age at recruitment/consent in your biobank
  
