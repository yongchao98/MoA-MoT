def analyze_aps_case():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the 2006 revised Sapporo criteria.
    """
    # --- Patient Data ---
    clinical_data = {
        "vte_events": 3,
        "notes": "Calf DVT at 18, PE at 25 (pregnancy), PE 4 months ago (spontaneous)"
    }
    
    # Lab results from 3 months ago and today
    lab_3_months_ago = {
        "antiB2GP1_IgM": 41, "antiB2GP1_IgG": 18,
        "anticardiolipin_IgM": 32, "anticardiolipin_IgG": 9,
        "dRVVT_ratio": 1.44
    }
    lab_today = {
        "antiB2GP1_IgM": 29, "antiB2GP1_IgG": 21,
        "anticardiolipin_IgM": 47, "anticardiolipin_IgG": 7,
        "dRVVT_ratio": 1.51
    }

    # --- Diagnostic Criteria Thresholds ---
    normal_values = {
        "antibody_cutoff": 20,      # UI/L, based on provided N < 20
        "dRVVT_ratio_cutoff": 1.2   # Based on provided N < 1.2
    }
    
    print("Evaluating for Antiphospholipid Syndrome (APS) based on revised Sapporo criteria.")
    print("Diagnosis requires at least one clinical and one laboratory criterion.\n")

    # --- Step 1: Assess Clinical Criteria ---
    print("--- Step 1: Assessing Clinical Criteria ---")
    clinical_criterion_met = False
    if clinical_data["vte_events"] >= 1:
        clinical_criterion_met = True
        print(f"Clinical Criterion: MET")
        print(f"Reason: Patient has a history of {clinical_data['vte_events']} documented venous thromboembolism (VTE) events.")
    else:
        print("Clinical Criterion: NOT MET")

    # --- Step 2: Assess Laboratory Criteria ---
    print("\n--- Step 2: Assessing Laboratory Criteria ---")
    print("A laboratory criterion must be confirmed on two occasions at least 12 weeks apart.")
    print("The patient's tests were performed 3 months ago and today, satisfying the interval.\n")

    # 2a. Lupus Anticoagulant (LA) based on dRVVT
    la_test1_positive = lab_3_months_ago['dRVVT_ratio'] > normal_values['dRVVT_ratio_cutoff']
    la_test2_positive = lab_today['dRVVT_ratio'] > normal_values['dRVVT_ratio_cutoff']
    persistent_la = la_test1_positive and la_test2_positive
    print("1. Lupus Anticoagulant (LA) Status:")
    print(f"  - 3 months ago (dRVVT ratio): {lab_3_months_ago['dRVVT_ratio']} (Normal < {normal_values['dRVVT_ratio_cutoff']}) -> Result: {'Positive' if la_test1_positive else 'Negative'}")
    print(f"  - Today (dRVVT ratio):        {lab_today['dRVVT_ratio']} (Normal < {normal_values['dRVVT_ratio_cutoff']}) -> Result: {'Positive' if la_test2_positive else 'Negative'}")
    if persistent_la:
        print("  -> Verdict: Persistent LA positivity confirmed. Criterion MET.")
    else:
        print("  -> Verdict: LA positivity not persistent. Criterion NOT MET.")
    print("     (Note: Rivaroxaban can cause false positive LA tests, but these consistently high values are strongly suggestive.)\n")

    # 2b. Anticardiolipin (aCL) Antibodies
    acl_test1_positive = (lab_3_months_ago['anticardiolipin_IgM'] > normal_values['antibody_cutoff'] or 
                          lab_3_months_ago['anticardiolipin_IgG'] > normal_values['antibody_cutoff'])
    acl_test2_positive = (lab_today['anticardiolipin_IgM'] > normal_values['antibody_cutoff'] or
                          lab_today['anticardiolipin_IgG'] > normal_values['antibody_cutoff'])
    persistent_acl = acl_test1_positive and acl_test2_positive
    print("2. Anticardiolipin (aCL) Antibody Status:")
    print(f"  - 3 months ago (IgM): {lab_3_months_ago['anticardiolipin_IgM']}, (IgG): {lab_3_months_ago['anticardiolipin_IgG']} (Normal < {normal_values['antibody_cutoff']}) -> Result: {'Positive' if acl_test1_positive else 'Negative'}")
    print(f"  - Today (IgM):        {lab_today['anticardiolipin_IgM']}, (IgG): {lab_today['anticardiolipin_IgG']} (Normal < {normal_values['antibody_cutoff']}) -> Result: {'Positive' if acl_test2_positive else 'Negative'}")
    if persistent_acl:
        print("  -> Verdict: Persistent aCL positivity confirmed. Criterion MET.\n")
    else:
        print("  -> Verdict: aCL positivity not persistent. Criterion NOT MET.\n")

    # 2c. Anti-β2 Glycoprotein-I (anti-β2GPI) Antibodies
    b2gp1_test1_positive = (lab_3_months_ago['antiB2GP1_IgM'] > normal_values['antibody_cutoff'] or 
                            lab_3_months_ago['antiB2GP1_IgG'] > normal_values['antibody_cutoff'])
    b2gp1_test2_positive = (lab_today['antiB2GP1_IgM'] > normal_values['antibody_cutoff'] or
                            lab_today['antiB2GP1_IgG'] > normal_values['antibody_cutoff'])
    persistent_b2gp1 = b2gp1_test1_positive and b2gp1_test2_positive
    print("3. Anti-β2 Glycoprotein-I (anti-β2GPI) Antibody Status:")
    print(f"  - 3 months ago (IgM): {lab_3_months_ago['antiB2GP1_IgM']}, (IgG): {lab_3_months_ago['antiB2GP1_IgG']} (Normal < {normal_values['antibody_cutoff']}) -> Result: {'Positive' if b2gp1_test1_positive else 'Negative'}")
    print(f"  - Today (IgM):        {lab_today['antiB2GP1_IgM']}, (IgG): {lab_today['antiB2GP1_IgG']} (Normal < {normal_values['antibody_cutoff']}) -> Result: {'Positive' if b2gp1_test2_positive else 'Negative'}")
    if persistent_b2gp1:
        print("  -> Verdict: Persistent anti-β2GPI positivity confirmed. Criterion MET.")
    else:
        print("  -> Verdict: anti-β2GPI positivity not persistent. Criterion NOT MET.")

    # --- Final Conclusion ---
    print("\n--- Final Diagnosis ---")
    laboratory_criterion_met = persistent_la or persistent_acl or persistent_b2gp1
    if clinical_criterion_met and laboratory_criterion_met:
        final_answer = "Yes"
        print("Conclusion: The patient meets both clinical (Vascular Thrombosis) and laboratory (persistently positive antiphospholipid antibodies) criteria.")
    else:
        final_answer = "No"
        print("Conclusion: The patient does not meet the full criteria for a definitive diagnosis of APS.")
    
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    analyze_aps_case()