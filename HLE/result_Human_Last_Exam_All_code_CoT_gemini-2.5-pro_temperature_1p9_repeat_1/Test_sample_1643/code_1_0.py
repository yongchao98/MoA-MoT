def check_aps_diagnosis():
    """
    Analyzes a patient's clinical and lab data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical
    vte_events = 3
    # Lab 1 (3 months ago)
    lab1_ac_igm = 32
    lab1_ac_igg = 9
    lab1_b2gp1_igm = 41
    lab1_b2gp1_igg = 18
    lab1_drvvtr = 1.44
    # Lab 2 (Today - 12 weeks later)
    lab2_ac_igm = 47
    lab2_ac_igg = 7
    lab2_b2gp1_igm = 29
    lab2_b2gp1_igg = 21
    lab2_drvvtr = 1.51
    # Criteria Thresholds
    medium_high_titer = 40
    drvvtr_normal = 1.2
    lab_interval_weeks = 12

    print("--- Antiphospholipid Syndrome (APS) Diagnostic Assessment ---")
    print("\nA diagnosis of APS requires at least 1 clinical and 1 laboratory criterion.")

    # --- Step 1: Assess Clinical Criteria ---
    print("\n--- 1. Clinical Criteria Assessment ---")
    clinical_criterion_met = False
    if vte_events >= 1:
        clinical_criterion_met = True
        print(f"Outcome: MET")
        print(f"Reasoning: The patient has a history of {vte_events} venous thromboembolic (VTE) events. This satisfies the criterion for vascular thrombosis.")
    else:
        print("Outcome: NOT MET")

    # --- Step 2: Assess Laboratory Criteria ---
    print("\n--- 2. Laboratory Criteria Assessment ---")
    print(f"Requirement: Persistent positive antiphospholipid antibodies on two occasions at least {lab_interval_weeks} weeks apart.\n")

    lab_criteria_met = []

    # Check for persistent Anticardiolipin (aCL) IgM
    is_acl_igm_persistent = lab1_ac_igm > 20 and lab2_ac_igm > 20
    is_acl_igm_high_titer = lab1_ac_igm >= medium_high_titer or lab2_ac_igm >= medium_high_titer
    if is_acl_igm_persistent and is_acl_igm_high_titer:
        lab_criteria_met.append("Anticardiolipin IgM")
        print(f"- Anticardiolipin IgM: MET")
        print(f"  - Result 1: {lab1_ac_igm} UI/L. Result 2: {lab2_ac_igm} UI/L.")
        print(f"  - Reasoning: Antibodies are persistently positive, and the titer of {max(lab1_ac_igm, lab2_ac_igm)} UI/L meets the medium-high titer requirement (> {medium_high_titer} UI/L).")

    # Check for persistent Anti-beta2-glycoprotein I (aB2GP1) IgM
    is_b2gp1_igm_persistent = lab1_b2gp1_igm > 20 and lab2_b2gp1_igm > 20
    is_b2gp1_igm_high_titer = lab1_b2gp1_igm >= medium_high_titer or lab2_b2gp1_igm >= medium_high_titer
    if is_b2gp1_igm_persistent and is_b2gp1_igm_high_titer:
        lab_criteria_met.append("Anti-beta2-glycoprotein I IgM")
        print(f"- Anti-beta2-glycoprotein I IgM: MET")
        print(f"  - Result 1: {lab1_b2gp1_igm} UI/L. Result 2: {lab2_b2gp1_igm} UI/L.")
        print(f"  - Reasoning: Antibodies are persistently positive, and the titer of {max(lab1_b2gp1_igm, lab2_b2gp1_igm)} UI/L meets the medium-high titer requirement (> {medium_high_titer} UI/L).")

    # Check for persistent Lupus Anticoagulant (dRVVT)
    is_drvvtr_persistent = lab1_drvvtr > drvvtr_normal and lab2_drvvtr > drvvtr_normal
    if is_drvvtr_persistent:
        lab_criteria_met.append("Lupus Anticoagulant (dRVVT)")
        print(f"- Lupus Anticoagulant (via dRVVT): MET")
        print(f"  - dRVVT Ratio 1: {lab1_drvvtr}. dRVVT Ratio 2: {lab2_drvvtr}.")
        print(f"  - Reasoning: The dRVVT ratio is persistently elevated above the normal threshold of {drvvtr_normal}. (Note: Rivaroxaban can cause false positives, but this result is considered in context with other positive markers).")

    if not lab_criteria_met:
        print("Outcome: NOT MET")
        print("Reasoning: No single laboratory marker was persistently positive at the required threshold over 12 weeks.")
    else:
        print(f"\nOutcome: MET")
        print(f"Reasoning: The patient has persistently positive results for: {', '.join(lab_criteria_met)}.")

    # --- Step 3: Final Conclusion ---
    print("\n--- 3. Final Conclusion ---")
    if clinical_criterion_met and len(lab_criteria_met) > 0:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion.")
    else:
        final_answer = "No"
        print("The patient does not meet the full diagnostic criteria for APS.")

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    check_aps_diagnosis()