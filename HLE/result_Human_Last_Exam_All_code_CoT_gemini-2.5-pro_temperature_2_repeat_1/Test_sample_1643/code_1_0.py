def analyze_aps_criteria():
    """
    Analyzes patient data to determine if they meet the criteria for Antiphospholipid Syndrome (APS)
    based on the 2006 revised Sapporo (Sydney) criteria.
    """

    # --- Step 1: Evaluate Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria")
    print("-----------------------------------")
    # The primary clinical criterion is "Vascular Thrombosis": one or more clinical episodes of
    # arterial, venous, or small vessel thrombosis.
    vte_events = 3
    clinical_criteria_met = vte_events >= 1
    
    if clinical_criteria_met:
        print(f"The patient has a history of {vte_events} venous thromboembolic (VTE) events.")
        print("Finding: The clinical criterion for 'Vascular Thrombosis' is MET.")
    else:
        print(f"The patient has {vte_events} documented VTE events.")
        print("Finding: The clinical criterion for 'Vascular Thrombosis' is NOT MET.")


    # --- Step 2: Evaluate Laboratory Criteria ---
    print("\nStep 2: Evaluating Laboratory Criteria")
    print("--------------------------------------")
    print("Criteria requires positivity on 2+ occasions at least 12 weeks apart.")

    # Lab data points
    lab_set_1 = {"dRVVT_ratio": 1.44, "aCL_IgM": 32, "aB2GPI_IgM": 41}
    lab_set_2 = {"dRVVT_ratio": 1.51, "aCL_IgM": 47, "aB2GPI_IgM": 29}
    
    # Normal/cutoff values
    normals = {"dRVVT_ratio_norm": 1.2, "aCL_IgM_medium_titer": 40, "aB2GPI_IgM_norm": 20}

    # 2a: Check for persistent Lupus Anticoagulant (LA) via dRVVT
    print("\n- Check 1: Lupus Anticoagulant (LA)")
    la_positive_1 = lab_set_1["dRVVT_ratio"] > normals["dRVVT_ratio_norm"]
    la_positive_2 = lab_set_2["dRVVT_ratio"] > normals["dRVVT_ratio_norm"]
    persistent_la = la_positive_1 and la_positive_2
    print(f"  - First test (dRVVT ratio): {lab_set_1['dRVVT_ratio']} (> {normals['dRVVT_ratio_norm']} is positive) -> Positive: {la_positive_1}")
    print(f"  - Second test (dRVVT ratio): {lab_set_2['dRVVT_ratio']} (> {normals['dRVVT_ratio_norm']} is positive) -> Positive: {la_positive_2}")
    if persistent_la:
        print("  - Finding: LA is persistently positive. This laboratory criterion is MET.")
    else:
        print("  - Finding: LA is NOT persistently positive.")

    # 2b: Check for persistent Anticardiolipin (aCL) antibody
    print("\n- Check 2: Anticardiolipin (aCL) antibody")
    # According to strict criteria, a medium-high titer (>40) must be present on both occasions.
    acl_positive_1 = lab_set_1["aCL_IgM"] > normals["aCL_IgM_medium_titer"]
    acl_positive_2 = lab_set_2["aCL_IgM"] > normals["aCL_IgM_medium_titer"]
    persistent_acl = acl_positive_1 and acl_positive_2
    print(f"  - First test (aCL IgM): {lab_set_1['aCL_IgM']} UI/L (> {normals['aCL_IgM_medium_titer']} for medium titer) -> Medium Titer: {acl_positive_1}")
    print(f"  - Second test (aCL IgM): {lab_set_2['aCL_IgM']} UI/L (> {normals['aCL_IgM_medium_titer']} for medium titer) -> Medium Titer: {acl_positive_2}")
    if persistent_acl:
         print("  - Finding: aCL IgM is persistently positive at medium titer.")
    else:
        print("  - Finding: aCL IgM is NOT persistently positive at medium titer by strict criteria.")
    
    # 2c: Check for persistent Anti-β2-glycoprotein-I (aβ2GPI) antibody
    print("\n- Check 3: Anti-β2-glycoprotein-I (aβ2GPI) antibody")
    ab2gpi_positive_1 = lab_set_1["aB2GPI_IgM"] > normals["aB2GPI_IgM_norm"]
    ab2gpi_positive_2 = lab_set_2["aB2GPI_IgM"] > normals["aB2GPI_IgM_norm"]
    persistent_ab2gpi = ab2gpi_positive_1 and ab2gpi_positive_2
    print(f"  - First test (aβ2GPI IgM): {lab_set_1['aB2GPI_IgM']} UI/L (> {normals['aB2GPI_IgM_norm']} is positive) -> Positive: {ab2gpi_positive_1}")
    print(f"  - Second test (aβ2GPI IgM): {lab_set_2['aB2GPI_IgM']} UI/L (> {normals['aB2GPI_IgM_norm']} is positive) -> Positive: {ab2gpi_positive_2}")
    if persistent_ab2gpi:
        print("  - Finding: aβ2GPI IgM is persistently positive. This laboratory criterion is MET.")
    else:
        print("  - Finding: aβ2GPI IgM is NOT persistently positive.")
        
    laboratory_criteria_met = persistent_la or persistent_acl or persistent_ab2gpi

    # --- Step 3: Final Conclusion ---
    print("\nStep 3: Final Conclusion")
    print("------------------------")
    print("Diagnosis requires at least ONE clinical AND at least ONE laboratory criterion to be met.")
    print(f"- Clinical Criteria Status: {'MET' if clinical_criteria_met else 'NOT MET'}")
    print(f"- Laboratory Criteria Status: {'MET' if laboratory_criteria_met else 'NOT MET'}")

    if clinical_criteria_met and laboratory_criteria_met:
        final_answer = "Yes"
        print("\nConclusion: The patient meets both clinical and laboratory criteria.")
    else:
        final_answer = "No"
        print("\nConclusion: The patient does not meet the full criteria for APS.")
        
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

analyze_aps_criteria()