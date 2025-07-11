def check_aps_diagnosis():
    """
    Analyzes patient data against APS diagnostic criteria and prints the evaluation.
    """

    # --- Patient Data ---
    clinical_events = {
        "VTE_events": [
            {"type": "calf DVT", "age": 18},
            {"type": "PE", "age": 25},
            {"type": "PE", "age": 34}
        ]
    }
    
    # Lab data from 3 months ago (Time 1) and today (Time 2)
    # The interval is ~12 weeks, which meets the criteria for re-testing.
    labs_time1 = {
        "antiB2GP1_IgM": 41, "antiB2GP1_IgG": 18,
        "anticardiolipin_IgM": 32, "anticardiolipin_IgG": 9,
        "dRVVT_ratio": 1.44
    }
    
    labs_time2 = {
        "antiB2GP1_IgM": 29, "antiB2GP1_IgG": 21,
        "anticardiolipin_IgM": 47, "anticardiolipin_IgG": 7,
        "dRVVT_ratio": 1.51
    }

    # --- Normal Ranges ---
    normal_ranges = {
        "antibody_level": 20,
        "dRVVT_ratio": 1.2
    }

    print("Evaluating patient for Antiphospholipid Syndrome (APS) based on the 2006 Sydney revised classification criteria.\n")

    # --- Step 1: Evaluate Clinical Criteria ---
    print("--- Step 1: Evaluating Clinical Criteria ---")
    num_vte_events = len(clinical_events["VTE_events"])
    clinical_criterion_met = num_vte_events >= 1
    
    print(f"The patient has a history of {num_vte_events} venous thromboembolic (VTE) events.")
    print("The presence of one or more episodes of venous, arterial, or small vessel thrombosis is a clinical criterion for APS.")
    if clinical_criterion_met:
        print("Result: Clinical criterion is MET.\n")
    else:
        print("Result: Clinical criterion is NOT MET.\n")

    # --- Step 2: Evaluate Laboratory Criteria ---
    print("--- Step 2: Evaluating Laboratory Criteria ---")
    print("Lab tests must be positive on two occasions at least 12 weeks apart. The tests were performed 3 months apart, meeting this requirement.\n")

    # Check for persistent Lupus Anticoagulant (via dRVVT)
    la_t1 = labs_time1["dRVVT_ratio"] > normal_ranges["dRVVT_ratio"]
    la_t2 = labs_time2["dRVVT_ratio"] > normal_ranges["dRVVT_ratio"]
    persistent_la = la_t1 and la_t2
    print(f"1. Lupus Anticoagulant (dRVVT):")
    print(f"   - Test 1 (3 months ago): {labs_time1['dRVVT_ratio']} (Normal < {normal_ranges['dRVVT_ratio']}) -> {'Positive' if la_t1 else 'Negative'}")
    print(f"   - Test 2 (Today): {labs_time2['dRVVT_ratio']} (Normal < {normal_ranges['dRVVT_ratio']}) -> {'Positive' if la_t2 else 'Negative'}")
    print(f"   - Persistently Positive: {'Yes' if persistent_la else 'No'}\n")

    # Check for persistent Anticardiolipin antibodies
    acl_igm_t1 = labs_time1["anticardiolipin_IgM"] > normal_ranges["antibody_level"]
    acl_igm_t2 = labs_time2["anticardiolipin_IgM"] > normal_ranges["antibody_level"]
    persistent_acl_igm = acl_igm_t1 and acl_igm_t2
    print(f"2. Anticardiolipin (aCL) IgM:")
    print(f"   - Test 1 (3 months ago): {labs_time1['anticardiolipin_IgM']} UI/L (Normal < {normal_ranges['antibody_level']}) -> {'Positive' if acl_igm_t1 else 'Negative'}")
    print(f"   - Test 2 (Today): {labs_time2['anticardiolipin_IgM']} UI/L (Normal < {normal_ranges['antibody_level']}) -> {'Positive' if acl_igm_t2 else 'Negative'}")
    print(f"   - Persistently Positive: {'Yes' if persistent_acl_igm else 'No'}\n")

    # Check for persistent Anti-ß2GP1 antibodies
    b2gp1_igm_t1 = labs_time1["antiB2GP1_IgM"] > normal_ranges["antibody_level"]
    b2gp1_igm_t2 = labs_time2["antiB2GP1_IgM"] > normal_ranges["antibody_level"]
    persistent_b2gp1_igm = b2gp1_igm_t1 and b2gp1_igm_t2
    print(f"3. Anti-ß2GP1 IgM:")
    print(f"   - Test 1 (3 months ago): {labs_time1['antiB2GP1_IgM']} UI/L (Normal < {normal_ranges['antibody_level']}) -> {'Positive' if b2gp1_igm_t1 else 'Negative'}")
    print(f"   - Test 2 (Today): {labs_time2['antiB2GP1_IgM']} UI/L (Normal < {normal_ranges['antibody_level']}) -> {'Positive' if b2gp1_igm_t2 else 'Negative'}")
    print(f"   - Persistently Positive: {'Yes' if persistent_b2gp1_igm else 'No'}\n")
    
    lab_criterion_met = persistent_la or persistent_acl_igm or persistent_b2gp1_igm
    
    if lab_criterion_met:
        print("Result: At least one laboratory criterion is MET.\n")
    else:
        print("Result: Laboratory criteria are NOT MET.\n")

    # --- Step 3: Final Diagnosis ---
    print("--- Step 3: Final Diagnosis ---")
    aps_diagnosed = clinical_criterion_met and lab_criterion_met
    
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {lab_criterion_met}")

    if aps_diagnosed:
        print("Conclusion: The patient meets both clinical and laboratory criteria for a definitive diagnosis of Antiphospholipid Syndrome.\n")
        final_answer = "Yes"
    else:
        print("Conclusion: The patient does not meet the full criteria for a definitive diagnosis of Antiphospholipid Syndrome.\n")
        final_answer = "No"

    print("Does this patient categorizes as having antiphospholipid syndrome ?")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    check_aps_diagnosis()