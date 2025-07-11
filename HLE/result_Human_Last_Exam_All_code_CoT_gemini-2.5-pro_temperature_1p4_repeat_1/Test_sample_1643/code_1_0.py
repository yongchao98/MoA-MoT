def diagnose_aps():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the revised Sapporo (Sydney) criteria.
    """

    # --- Patient Data ---
    # Clinical
    vte_events = 3

    # Lab tests 3 months ago
    lab1 = {
        'anti_b2gp1_igm': 41, 'anti_b2gp1_igg': 18,
        'anticardiolipin_igm': 32, 'anticardiolipin_igg': 9,
        'drvvT_ratio': 1.44
    }

    # Lab tests today (3 months later, satisfying the >= 12 weeks interval)
    lab2 = {
        'anti_b2gp1_igm': 29, 'anti_b2gp1_igg': 21,
        'anticardiolipin_igm': 47, 'anticardiolipin_igg': 7,
        'drvvT_ratio': 1.51
    }

    # Normal ranges / Positive thresholds
    normal_antibody_upper = 20  # UI/L
    normal_drvvT_ratio_upper = 1.2

    print("--- Antiphospholipid Syndrome (APS) Diagnostic Check ---")

    # --- Step 1: Evaluate Clinical Criteria ---
    # Criterion: One or more clinical episodes of arterial, venous, or small vessel thrombosis.
    clinical_criterion_met = False
    print("\n1. Evaluating Clinical Criteria:")
    if vte_events >= 1:
        clinical_criterion_met = True
        print(f"   - Patient has a history of {vte_events} venous thromboembolism (VTE) events.")
        print("   - Result: Clinical criterion is MET.")
    else:
        print("   - Patient does not have a documented history of thrombosis.")
        print("   - Result: Clinical criterion is NOT MET.")

    # --- Step 2: Evaluate Laboratory Criteria ---
    # Criterion: At least one positive lab test, confirmed on two occasions >= 12 weeks apart.
    lab_criterion_met = False
    la_positive = False
    acl_positive = False
    ab2gp1_positive = False

    print("\n2. Evaluating Laboratory Criteria (must be persistent over >= 12 weeks):")

    # Check for Lupus Anticoagulant (LA) via dRVVT
    # Note: Rivaroxaban (a DOAC) can cause false positives in dRVVT tests.
    # However, for this evaluation, we will assess the given data.
    print("\n   a) Lupus Anticoagulan (via dRVVT):")
    if lab1['drvvT_ratio'] > normal_drvvT_ratio_upper and lab2['drvvT_ratio'] > normal_drvvT_ratio_upper:
        la_positive = True
        print(f"      - First dRVVT ratio was {lab1['drvvT_ratio']} (Normal < {normal_drvvT_ratio_upper}).")
        print(f"      - Second dRVVT ratio was {lab2['drvvT_ratio']} (Normal < {normal_drvvT_ratio_upper}).")
        print("      - Result: Persistently positive. This sub-criterion is MET.")
        print("      - Clinical Note: Patient is on Rivaroxaban, which can cause false positives for this test.")
    else:
        print("      - Result: Not persistently positive. This sub-criterion is NOT MET.")

    # Check for anticardiolipin (aCL) antibodies
    # Criterion is medium-to-high titer (>40 GPL/MPL or >99th percentile).
    # We will consider any persistent positive result as significant here, especially with the 2nd value being >40.
    print("\n   b) Anticardiolipin Antibodies:")
    if (lab1['anticardiolipin_igm'] > normal_antibody_upper and lab2['anticardiolipin_igm'] > normal_antibody_upper) or \
       (lab1['anticardiolipin_igg'] > normal_antibody_upper and lab2['anticardiolipin_igg'] > normal_antibody_upper):
        acl_positive = True
        print(f"      - Anticardiolipin IgM was {lab1['anticardiolipin_igm']} then {lab2['anticardiolipin_igm']} (Normal < {normal_antibody_upper}).")
        print(f"      - Anticardiolipin IgG was {lab1['anticardiolipin_igg']} then {lab2['anticardiolipin_igg']} (Normal < {normal_antibody_upper}).")
        print("      - Result: Persistently positive aCL IgM. This sub-criterion is MET.")
    else:
        print("      - Result: Not persistently positive. This sub-criterion is NOT MET.")

    # Check for anti-beta-2-glycoprotein-I (aB2GPI) antibodies
    print("\n   c) Anti-ß2-Glycoprotein-I Antibodies:")
    if (lab1['anti_b2gp1_igm'] > normal_antibody_upper and lab2['anti_b2gp1_igm'] > normal_antibody_upper) or \
       (lab1['anti_b2gp1_igg'] > normal_antibody_upper and lab2['anti_b2gp1_igg'] > normal_antibody_upper):
        ab2gp1_positive = True
        print(f"      - Anti-ß2GP1 IgM was {lab1['anti_b2gp1_igm']} then {lab2['anti_b2gp1_igm']} (Normal < {normal_antibody_upper}).")
        # Check for newly positive IgG
        if lab1['anti_b2gp1_igg'] <= normal_antibody_upper and lab2['anti_b2gp1_igg'] > normal_antibody_upper:
             print(f"      - Anti-ß2GP1 IgG became positive: {lab1['anti_b2gp1_igg']} then {lab2['anti_b2gp1_igg']} (Normal < {normal_antibody_upper}).")
        print("      - Result: Persistently positive anti-ß2GP1 IgM. This sub-criterion is MET.")
    else:
        print("      - Result: Not persistently positive. This sub-criterion is NOT MET.")

    if la_positive or acl_positive or ab2gp1_positive:
        lab_criterion_met = True
        print("\n   - Overall Lab Result: At least one laboratory criterion is MET.")
    else:
        print("\n   - Overall Lab Result: No laboratory criteria are met.")


    # --- Step 3: Final Conclusion ---
    print("\n--- CONCLUSION ---")
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets both clinical AND laboratory criteria for Antiphospholipid Syndrome.")
    else:
        final_answer = "No"
        print("The patient does not meet the full criteria for Antiphospholipid Syndrome.")
    
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

diagnose_aps()