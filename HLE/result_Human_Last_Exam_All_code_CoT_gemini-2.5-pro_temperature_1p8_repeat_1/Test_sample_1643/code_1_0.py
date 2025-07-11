def check_aps_diagnosis():
    """
    Evaluates a patient's clinical and laboratory data against the Sydney criteria
    for Antiphospholipid Syndrome (APS) and prints the reasoning.
    """

    # --- Patient Data ---
    patient_info = {
        "vte_events": 3,
    }

    lab_results = {
        "test_1": {"interval_months": 3, "dRVVT_ratio": 1.44, "aCL_IgM": 32, "aB2GP1_IgM": 41, "aB2GP1_IgG": 18},
        "test_2": {"interval_months": 0, "dRVVT_ratio": 1.51, "aCL_IgM": 47, "aB2GP1_IgM": 29, "aB2GP1_IgG": 21},
    }

    # --- APS Diagnostic Criteria ---
    criteria = {
        "clinical_thrombosis_events": 1,
        "lab_time_apart_weeks": 12,
        "dRVVT_ratio_cutoff": 1.2,
        "aCL_cutoff": 20,
        "aCL_medium_titer_cutoff": 40,
        "aB2GP1_cutoff": 20,
    }

    print("--- Evaluating for Antiphospholipid Syndrome (APS) ---\n")

    # --- 1. Clinical Criteria Evaluation ---
    clinical_met = False
    print("Step 1: Evaluating Clinical Criteria (requires >=1 vascular thrombosis event)")
    if patient_info["vte_events"] >= criteria["clinical_thrombosis_events"]:
        clinical_met = True
        print(f"Result: Clinical criterion MET.")
        print(f"Reason: The patient has had {patient_info['vte_events']} VTE events, which meets the requirement of >= {criteria['clinical_thrombosis_events']} event.\n")
    else:
        print(f"Result: Clinical criterion NOT MET.\n")

    # --- 2. Laboratory Criteria Evaluation ---
    lab_met = False
    print("Step 2: Evaluating Laboratory Criteria (requires >=1 positive marker, confirmed >12 weeks apart)")

    # A. Lupus Anticoagulant (LA) based on dRVVT
    lab1 = lab_results["test_1"]
    lab2 = lab_results["test_2"]
    
    la_positive_1 = lab1["dRVVT_ratio"] > criteria["dRVVT_ratio_cutoff"]
    la_positive_2 = lab2["dRVVT_ratio"] > criteria["dRVVT_ratio_cutoff"]
    la_persistent = la_positive_1 and la_positive_2

    print(f"- Checking for persistent Lupus Anticoagulant...")
    print(f"  Test 1 dRVVT ratio: {lab1['dRVVT_ratio']} (Normal < {criteria['dRVVT_ratio_cutoff']}) -> {'Positive' if la_positive_1 else 'Negative'}")
    print(f"  Test 2 dRVVT ratio: {lab2['dRVVT_ratio']} (Normal < {criteria['dRVVT_ratio_cutoff']}) -> {'Positive' if la_positive_2 else 'Negative'}")
    if la_persistent:
        lab_met = True
        print(f"  Result: LA criterion MET.\n")
    else:
        print(f"  Result: LA criterion NOT MET.\n")

    # B. Anticardiolipin (aCL) antibodies
    acl_igm_positive_1 = lab1["aCL_IgM"] > criteria["aCL_cutoff"]
    acl_igm_positive_2 = lab2["aCL_IgM"] > criteria["aCL_cutoff"]
    # Medium titer is > 40
    acl_medium_titer_reached = lab1["aCL_IgM"] > criteria["aCL_medium_titer_cutoff"] or lab2["aCL_IgM"] > criteria["aCL_medium_titer_cutoff"]
    acl_persistent_and_medium = acl_igm_positive_1 and acl_igm_positive_2 and acl_medium_titer_reached
    
    print(f"- Checking for persistent Anticardiolipin (aCL) antibodies at medium/high titer...")
    print(f"  Test 1 aCL IgM: {lab1['aCL_IgM']} UI/L (Normal < {criteria['aCL_cutoff']}) -> {'Positive' if acl_igm_positive_1 else 'Negative'}")
    print(f"  Test 2 aCL IgM: {lab2['aCL_IgM']} UI/L (Normal < {criteria['aCL_cutoff']}) -> {'Positive' if acl_igm_positive_2 else 'Negative'}")
    if acl_persistent_and_medium:
        lab_met = True # This would also satisfy the lab criteria
        print(f"  Result: aCL criterion MET (persistently positive and reached medium titer of >{criteria['aCL_medium_titer_cutoff']}).\n")
    else:
        print(f"  Result: aCL criterion NOT MET.\n")


    # C. Anti-Beta2-Glycoprotein-I (aB2GP1) antibodies
    ab2gp1_igm_positive_1 = lab1['aB2GP1_IgM'] > criteria['aB2GP1_cutoff']
    ab2gp1_igm_positive_2 = lab2['aB2GP1_IgM'] > criteria['aB2GP1_cutoff']
    ab2gp1_igm_persistent = ab2gp1_igm_positive_1 and ab2gp1_igm_positive_2
    
    print(f"- Checking for persistent Anti-Beta2-Glycoprotein-I (aB2GP1) antibodies...")
    print(f"  Test 1 aB2GP1 IgM: {lab1['aB2GP1_IgM']} UI/L (Normal < {criteria['aB2GP1_cutoff']}) -> {'Positive' if ab2gp1_igm_positive_1 else 'Negative'}")
    print(f"  Test 2 aB2GP1 IgM: {lab2['aB2GP1_IgM']} UI/L (Normal < {criteria['aB2GP1_cutoff']}) -> {'Positive' if ab2gp1_igm_positive_2 else 'Negative'}")
    if ab2gp1_igm_persistent:
        lab_met = True # This would also satisfy the lab criteria
        print(f"  Result: aB2GP1 criterion MET.\n")
    else:
        print(f"  Result: aB2GP1 criterion NOT MET.\n")
    
    # --- 3. Final Conclusion ---
    print("--- Conclusion ---")
    if clinical_met and lab_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion AND at least one laboratory criterion.")
        print("Diagnosis of Antiphospholipid Syndrome is confirmed.")
    else:
        final_answer = "No"
        print("The patient does not meet the full criteria for an APS diagnosis.")
    
    print("\nDoes this patient categorizes as having antiphospholipid syndrome ?")
    print(f"<<<{final_answer}>>>")

# Run the evaluation
check_aps_diagnosis()