def check_antiphospholipid_syndrome():
    """
    Evaluates a patient's clinical and lab data against the diagnostic
    criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical
    vte_events = 3

    # Lab - Time 1 (3 months ago)
    dRVVT_ratio_1 = 1.44
    acl_igm_1 = 32
    ab2gp1_igm_1 = 41

    # Lab - Time 2 (Today)
    dRVVT_ratio_2 = 1.51
    acl_igm_2 = 47
    ab2gp1_igm_2 = 29
    
    # Lab - Time 2 (Today) - other values for completeness
    ab2gp1_igg_2 = 21

    # --- APS Criteria Cutoffs ---
    dRVVT_cutoff = 1.2
    antibody_normal_cutoff = 20
    antibody_medium_titer_cutoff = 40

    print("Evaluating for Antiphospholipid Syndrome (APS)...")
    print("-" * 50)

    # 1. Check Clinical Criteria
    print("Step 1: Assessing Clinical Criteria")
    meets_clinical_criteria = vte_events >= 1
    print(f"Patient has {vte_events} VTE events. Required: >= 1.")
    if meets_clinical_criteria:
        print("Result: Clinical criteria are MET.\n")
    else:
        print("Result: Clinical criteria are NOT MET.\n")

    # 2. Check Laboratory Criteria
    print("Step 2: Assessing Laboratory Criteria (Persistence over >=12 weeks)")

    # Check for persistent Lupus Anticoagulant (LA) via dRVVT
    la_positive_t1 = dRVVT_ratio_1 > dRVVT_cutoff
    la_positive_t2 = dRVVT_ratio_2 > dRVVT_cutoff
    persistent_la = la_positive_t1 and la_positive_t2
    print(f"Lupus Anticoagulant (dRVVT > {dRVVT_cutoff}):")
    print(f"  - Test 1: {dRVVT_ratio_1} -> {'Positive' if la_positive_t1 else 'Negative'}")
    print(f"  - Test 2: {dRVVT_ratio_2} -> {'Positive' if la_positive_t2 else 'Negative'}")
    print(f"  - Persistent LA: {'Yes' if persistent_la else 'No'}")

    # Check for persistent Anticardiolipin (aCL) antibodies (IgM in this case)
    acl_positive_t1 = acl_igm_1 > antibody_normal_cutoff
    acl_positive_t2 = acl_igm_2 > antibody_normal_cutoff
    acl_medium_titer = acl_igm_1 > antibody_medium_titer_cutoff or acl_igm_2 > antibody_medium_titer_cutoff
    persistent_acl = (acl_positive_t1 and acl_positive_t2 and acl_medium_titer)
    print(f"\nAnticardiolipin IgM (> {antibody_normal_cutoff} U/L and medium titer > {antibody_medium_titer_cutoff} U/L):")
    print(f"  - Test 1: {acl_igm_1} -> {'Positive' if acl_positive_t1 else 'Negative'}")
    print(f"  - Test 2: {acl_igm_2} -> {'Positive' if acl_positive_t2 else 'Negative'}")
    print(f"  - Persistent aCL IgM with medium titer: {'Yes' if persistent_acl else 'No'}")
    
    # Check for persistent Anti-beta2-glycoprotein-I (aB2GPI) antibodies (IgM in this case)
    ab2gp1_positive_t1 = ab2gp1_igm_1 > antibody_normal_cutoff
    ab2gp1_positive_t2 = ab2gp1_igm_2 > antibody_normal_cutoff or ab2gp1_igg_2 > antibody_normal_cutoff
    ab2gp1_medium_titer = ab2gp1_igm_1 > antibody_medium_titer_cutoff
    persistent_ab2gp1 = (ab2gp1_positive_t1 and ab2gp1_positive_t2 and ab2gp1_medium_titer)
    print(f"\nAnti-B2GP1 IgM (> {antibody_normal_cutoff} U/L and medium titer > {antibody_medium_titer_cutoff} U/L):")
    print(f"  - Test 1: {ab2gp1_igm_1} -> {'Positive' if ab2gp1_positive_t1 else 'Negative'}")
    print(f"  - Test 2: {ab2gp1_igm_2} (IgM) or {ab2gp1_igg_2} (IgG) -> {'Positive' if ab2gp1_positive_t2 else 'Negative'}")
    print(f"  - Persistent aB2GP1 with medium titer: {'Yes' if persistent_ab2gp1 else 'No'}")

    meets_lab_criteria = persistent_la or persistent_acl or persistent_ab2gp1
    print("\nResult: At least one laboratory criterion is persistently positive.")
    if meets_lab_criteria:
        print("Overall laboratory criteria are MET.\n")
    else:
        print("Overall laboratory criteria are NOT MET.\n")

    # 3. Final Conclusion
    print("-" * 50)
    print("Step 3: Final Conclusion")
    is_aps = meets_clinical_criteria and meets_lab_criteria
    final_answer = "Yes" if is_aps else "No"
    
    print(f"Clinical criteria met: {meets_clinical_criteria}")
    print(f"Laboratory criteria met: {meets_lab_criteria}")
    print(f"Diagnosis of Antiphospholipid Syndrome: {final_answer}")
    
    return final_answer

if __name__ == "__main__":
    answer = check_antiphospholipid_syndrome()
    print(f"<<<{answer}>>>")