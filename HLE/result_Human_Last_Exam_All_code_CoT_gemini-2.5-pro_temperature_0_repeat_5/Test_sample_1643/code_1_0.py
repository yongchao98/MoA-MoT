import sys

# This script is for educational and illustrative purposes and does not constitute medical advice.
# The interference of anticoagulants like Rivaroxaban on coagulation-based assays (e.g., dRVVT)
# should be considered in a real clinical setting.

def solve():
    """
    Analyzes patient data against the Sydney classification criteria for Antiphospholipid Syndrome (APS).
    """
    # --- Patient Data ---

    # Clinical Data
    vte_events = 3  # DVT + 2 PEs

    # Lab Data - Timepoint 1 (3 months ago)
    t1_dRVVT_ratio = 1.44
    t1_dRVVT_norm = 1.2
    t1_ac_igm = 32
    t1_ac_igg = 9
    t1_b2gp1_igm = 41
    t1_b2gp1_igg = 18
    t1_ac_norm = 20
    t1_b2gp1_norm = 20

    # Lab Data - Timepoint 2 (Today)
    # Interval is 3 months, which is >= 12 weeks.
    t2_dRVVT_ratio = 1.51
    t2_dRVVT_norm = 1.2
    t2_ac_igm = 47
    t2_ac_igg = 7
    t2_b2gp1_igm = 29
    t2_b2gp1_igg = 21
    t2_ac_norm = 20
    t2_b2gp1_norm = 20

    # --- Analysis ---

    print("Analyzing patient data for Antiphospholipid Syndrome (APS) diagnosis...")
    print("-" * 30)

    # 1. Check Clinical Criteria
    print("Step 1: Evaluating Clinical Criteria (requires >= 1 vascular thrombosis event)")
    clinical_criterion_met = vte_events >= 1
    print(f"Result: Patient has {vte_events} VTE events. Criterion met: {clinical_criterion_met}")
    print("-" * 30)

    # 2. Check Laboratory Criteria
    print("Step 2: Evaluating Laboratory Criteria (requires >= 1 persistent positive antibody test > 12 weeks apart)")

    # 2a. Lupus Anticoagulant (LA) - using dRVVT as the primary indicator here
    la_positive_t1 = t1_dRVVT_ratio > t1_dRVVT_norm
    la_positive_t2 = t2_dRVVT_ratio > t2_dRVVT_norm
    persistent_la = la_positive_t1 and la_positive_t2
    print("  - Lupus Anticoagulant (LA) based on dRVVT:")
    print(f"    Test 1: dRVVT ratio is {t1_dRVVT_ratio} (Normal < {t1_dRVVT_norm}). Positive: {la_positive_t1}")
    print(f"    Test 2: dRVVT ratio is {t2_dRVVT_ratio} (Normal < {t2_dRVVT_norm}). Positive: {la_positive_t2}")
    print(f"    Criterion met (persistent LA): {persistent_la}")
    print()

    # 2b. Anticardiolipin (aCL) Antibodies
    # Must be persistent and medium/high titer (>40 IU/L or >99th percentile)
    acl_igm_positive_t1 = t1_ac_igm > t1_ac_norm
    acl_igm_positive_t2 = t2_ac_igm > t2_ac_norm
    acl_igg_positive_t1 = t1_ac_igg > t1_ac_norm
    acl_igg_positive_t2 = t2_ac_igg > t2_ac_norm
    
    persistent_acl_igm = acl_igm_positive_t1 and acl_igm_positive_t2
    persistent_acl_igg = acl_igg_positive_t1 and acl_igg_positive_t2
    
    medium_high_titer_present = (t1_ac_igm > 40) or (t2_ac_igm > 40) or (t1_ac_igg > 40) or (t2_ac_igg > 40)
    
    persistent_acl = (persistent_acl_igm or persistent_acl_igg) and medium_high_titer_present
    print("  - Anticardiolipin (aCL) Antibodies:")
    print(f"    Test 1: aCL IgM {t1_ac_igm}, aCL IgG {t1_ac_igg} (Normal < {t1_ac_norm})")
    print(f"    Test 2: aCL IgM {t2_ac_igm}, aCL IgG {t2_ac_igg} (Normal < {t2_ac_norm})")
    print(f"    Persistent positive aCL (IgM or IgG): {persistent_acl_igm or persistent_acl_igg}")
    print(f"    Medium/High Titer (>40) present: {medium_high_titer_present} (since {t2_ac_igm} > 40)")
    print(f"    Criterion met (persistent aCL at medium/high titer): {persistent_acl}")
    print()

    # 2c. Anti-β2-glycoprotein I (anti-β2GPI) Antibodies
    # Must be persistent and >99th percentile (using >20 as the threshold)
    b2gp1_igm_positive_t1 = t1_b2gp1_igm > t1_b2gp1_norm
    b2gp1_igm_positive_t2 = t2_b2gp1_igm > t2_b2gp1_norm
    b2gp1_igg_positive_t1 = t1_b2gp1_igg > t1_b2gp1_norm
    b2gp1_igg_positive_t2 = t2_b2gp1_igg > t2_b2gp1_norm

    persistent_b2gp1_igm = b2gp1_igm_positive_t1 and b2gp1_igm_positive_t2
    persistent_b2gp1_igg = b2gp1_igg_positive_t1 and b2gp1_igg_positive_t2

    persistent_b2gp1 = persistent_b2gp1_igm or persistent_b2gp1_igg
    print("  - Anti-β2-glycoprotein I (anti-β2GPI) Antibodies:")
    print(f"    Test 1: anti-β2GPI IgM {t1_b2gp1_igm}, anti-β2GPI IgG {t1_b2gp1_igg} (Normal < {t1_b2gp1_norm})")
    print(f"    Test 2: anti-β2GPI IgM {t2_b2gp1_igm}, anti-β2GPI IgG {t2_b2gp1_igg} (Normal < {t2_b2gp1_norm})")
    print(f"    Persistent positive anti-β2GPI (IgM or IgG): {persistent_b2gp1}")
    print(f"    Criterion met (persistent anti-β2GPI): {persistent_b2gp1}")
    print("-" * 30)

    # Overall Lab Criteria
    lab_criteria_met = persistent_la or persistent_acl or persistent_b2gp1
    print(f"Step 2 Result: At least one laboratory criterion was met: {lab_criteria_met}")
    print("-" * 30)

    # 3. Final Conclusion
    print("Step 3: Final Conclusion (requires both clinical and lab criteria to be met)")
    is_aps = clinical_criterion_met and lab_criteria_met
    
    final_answer = "No"
    if is_aps:
        final_answer = "Yes"
        
    print(f"Clinical criterion met: {clinical_criterion_met}")
    print(f"Laboratory criterion met: {lab_criteria_met}")
    print(f"Does the patient categorize as having antiphospholipid syndrome? {final_answer}")
    
    # Final answer in the required format
    print(f"<<<{final_answer}>>>")

solve()