import sys

def diagnose_aps():
    """
    Analyzes patient data against the Revised Sapporo (Sydney) criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical History
    vte_events = 3  # Calf DVT, PE, another PE

    # Lab tests 3 months ago
    lab_1 = {
        'anti_b2gp1_igm': 41, 'anti_b2gp1_igg': 18,
        'anticardiolipin_igm': 32,
        'drvvt_ratio': 1.44
    }

    # Lab tests today (12 weeks later)
    lab_2 = {
        'anti_b2gp1_igm': 29, 'anti_b2gp1_igg': 21,
        'anticardiolipin_igm': 47,
        'drvvt_ratio': 1.51
    }

    # --- APS Criteria Thresholds ---
    VTE_MIN_EVENTS = 1
    ANTIBODY_NORMAL_CUTOFF = 20
    ACL_MEDIUM_HIGH_TITER = 40
    DRVVT_NORMAL_CUTOFF = 1.2
    
    print("--- Evaluating Antiphospholipid Syndrome (APS) Criteria ---")
    
    # 1. Evaluate Clinical Criteria
    print("\nStep 1: Checking Clinical Criteria...")
    clinical_criterion_met = vte_events >= VTE_MIN_EVENTS
    print(f"Criterion: Vascular Thrombosis (>= {VTE_MIN_EVENTS} event).")
    print(f"Patient's VTE Events: {vte_events}. Is criterion met? {clinical_criterion_met}")

    # 2. Evaluate Laboratory Criteria (must be persistent over 12 weeks)
    print("\nStep 2: Checking for Persistent Laboratory Criteria...")
    
    # a) Lupus Anticoagulant (LA) via dRVVT
    la_positive_t1 = lab_1['drvvt_ratio'] > DRVVT_NORMAL_CUTOFF
    la_positive_t2 = lab_2['drvvt_ratio'] > DRVVT_NORMAL_CUTOFF
    persistent_la = la_positive_t1 and la_positive_t2
    print(f"\n- Criterion: Persistent Lupus Anticoagulant (dRVVT > {DRVVT_NORMAL_CUTOFF})")
    print(f"  dRVVT 3 months ago: {lab_1['drvvt_ratio']}. Positive? {la_positive_t1}")
    print(f"  dRVVT Today: {lab_2['drvvt_ratio']}. Positive? {la_positive_t2}")
    print(f"  Is criterion met? {persistent_la}")
    print("  (Note: dRVVT may be falsely elevated by Rivaroxaban treatment)")

    # b) Anticardiolipin (aCL) Antibodies
    acl_positive_t1 = lab_1['anticardiolipin_igm'] > ANTIBODY_NORMAL_CUTOFF
    acl_positive_t2 = lab_2['anticardiolipin_igm'] > ANTIBODY_NORMAL_CUTOFF
    # Criterion requires medium/high titer (>40) on at least one of the positive occasions
    acl_has_high_titer = lab_1['anticardiolipin_igm'] > ACL_MEDIUM_HIGH_TITER or \
                         lab_2['anticardiolipin_igm'] > ACL_MEDIUM_HIGH_TITER
    persistent_acl = acl_positive_t1 and acl_positive_t2 and acl_has_high_titer
    print(f"\n- Criterion: Persistent Anticardiolipin IgM Antibodies (> {ANTIBODY_NORMAL_CUTOFF} U/mL and titer > {ACL_MEDIUM_HIGH_TITER} on at least one occasion)")
    print(f"  aCL IgM 3 months ago: {lab_1['anticardiolipin_igm']}. Positive? {acl_positive_t1}")
    print(f"  aCL IgM Today: {lab_2['anticardiolipin_igm']}. Positive? {acl_positive_t2}. Titer > {ACL_MEDIUM_HIGH_TITER}? {lab_2['anticardiolipin_igm'] > ACL_MEDIUM_HIGH_TITER}")
    print(f"  Is criterion met? {persistent_acl}")

    # c) Anti-beta2-glycoprotein-I (aB2GP1) Antibodies
    ab2gp1_positive_t1 = lab_1['anti_b2gp1_igm'] > ANTIBODY_NORMAL_CUTOFF
    ab2gp1_positive_t2 = lab_2['anti_b2gp1_igm'] > ANTIBODY_NORMAL_CUTOFF or lab_2['anti_b2gp1_igg'] > ANTIBODY_NORMAL_CUTOFF
    persistent_ab2gp1 = ab2gp1_positive_t1 and ab2gp1_positive_t2
    print(f"\n- Criterion: Persistent Anti-ß2GP1 Antibodies (> {ANTIBODY_NORMAL_CUTOFF} U/mL)")
    print(f"  Anti-ß2GP1 IgM 3 months ago: {lab_1['anti_b2gp1_igm']}. Positive? {ab2gp1_positive_t1}")
    print(f"  Anti-ß2GP1 IgM/IgG Today: IgM {lab_2['anti_b2gp1_igm']}, IgG {lab_2['anti_b2gp1_igg']}. Positive? {ab2gp1_positive_t2}")
    print(f"  Is criterion met? {persistent_ab2gp1}")
    
    # At least one lab criterion must be met
    lab_criteria_met = persistent_la or persistent_acl or persistent_ab2gp1

    # 3. Final Diagnosis
    print("\n--- Final Conclusion ---")
    is_aps = clinical_criterion_met and lab_criteria_met
    
    print(f"Final Equation: Is APS = (Clinical Criterion Met) AND (At least one Lab Criterion Met)")
    print(f"Result: {clinical_criterion_met} AND {lab_criteria_met} = {is_aps}")
    
    if is_aps:
        print("\nDoes this patient categorize as having antiphospholipid syndrome? Yes.")
        # Final answer format for the system
        sys.stdout.write("<<<Yes>>>")
    else:
        print("\nDoes this patient categorize as having antiphospholipid syndrome? No.")
        # Final answer format for the system
        sys.stdout.write("<<<No>>>")

if __name__ == '__main__':
    diagnose_aps()