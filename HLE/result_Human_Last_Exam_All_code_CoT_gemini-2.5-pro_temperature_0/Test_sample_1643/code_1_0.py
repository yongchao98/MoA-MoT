def diagnose_aps():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the 2006 Sydney consensus criteria.
    """
    print("Evaluating patient for Antiphospholipid Syndrome (APS)...")
    print("-" * 50)

    # --- 1. Clinical Criteria Evaluation ---
    print("Step 1: Assessing Clinical Criteria")
    vte_event_count = 3
    print(f"Patient has a history of {vte_event_count} venous thromboembolic (VTE) events.")
    
    # APS requires at least one clinical episode of vascular thrombosis.
    clinical_criterion_met = vte_event_count >= 1
    if clinical_criterion_met:
        print("Result: Clinical criterion for vascular thrombosis is MET.\n")
    else:
        print("Result: Clinical criterion for vascular thrombosis is NOT MET.\n")

    # --- 2. Laboratory Criteria Evaluation ---
    print("Step 2: Assessing Laboratory Criteria (Persistence over >=12 weeks)")
    # Lab data from two time points, 3 months apart.
    
    # Lupus Anticoagulant (LA) via dRVVT
    drvvt_t1 = 1.44
    drvvt_t2 = 1.51
    drvvt_normal = 1.2
    persistent_la = drvvt_t1 > drvvt_normal and drvvt_t2 > drvvt_normal
    print(f"Lupus Anticoagulant (dRVVT):")
    print(f"  - Test 1 (3 months ago): {drvvt_t1} (Normal < {drvvt_normal}) -> {'Positive' if drvvt_t1 > drvvt_normal else 'Negative'}")
    print(f"  - Test 2 (Today): {drvvt_t2} (Normal < {drvvt_normal}) -> {'Positive' if drvvt_t2 > drvvt_normal else 'Negative'}")
    print(f"  - Persistence Check: {'MET' if persistent_la else 'NOT MET'}\n")

    # Anticardiolipin (aCL) IgM
    acl_igm_t1 = 32
    acl_igm_t2 = 47
    acl_normal = 20
    persistent_acl = acl_igm_t1 > acl_normal and acl_igm_t2 > acl_normal
    print(f"Anticardiolipin IgM:")
    print(f"  - Test 1 (3 months ago): {acl_igm_t1} UI/L (Normal < {acl_normal}) -> {'Positive' if acl_igm_t1 > acl_normal else 'Negative'}")
    print(f"  - Test 2 (Today): {acl_igm_t2} UI/L (Normal < {acl_normal}) -> {'Positive' if acl_igm_t2 > acl_normal else 'Negative'}")
    print(f"  - Persistence Check: {'MET' if persistent_acl else 'NOT MET'}\n")

    # Anti-beta2-glycoprotein-I (aB2GP1) IgM
    ab2gp1_igm_t1 = 41
    ab2gp1_igm_t2 = 29
    ab2gp1_normal = 20
    persistent_ab2gp1 = ab2gp1_igm_t1 > ab2gp1_normal and ab2gp1_igm_t2 > ab2gp1_normal
    print(f"Anti-B2GP1 IgM:")
    print(f"  - Test 1 (3 months ago): {ab2gp1_igm_t1} UI/L (Normal < {ab2gp1_normal}) -> {'Positive' if ab2gp1_igm_t1 > ab2gp1_normal else 'Negative'}")
    print(f"  - Test 2 (Today): {ab2gp1_igm_t2} UI/L (Normal < {ab2gp1_normal}) -> {'Positive' if ab2gp1_igm_t2 > ab2gp1_normal else 'Negative'}")
    print(f"  - Persistence Check: {'MET' if persistent_ab2gp1 else 'NOT MET'}\n")

    # APS requires at least one positive laboratory criterion.
    laboratory_criterion_met = persistent_la or persistent_acl or persistent_ab2gp1
    if laboratory_criterion_met:
        print("Result: At least one laboratory criterion is MET.\n")
    else:
        print("Result: Laboratory criteria are NOT MET.\n")

    # --- 3. Final Conclusion ---
    print("-" * 50)
    print("Step 3: Final Diagnosis")
    
    aps_diagnosis = clinical_criterion_met and laboratory_criterion_met
    
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {laboratory_criterion_met}")
    
    if aps_diagnosis:
        print("Conclusion: The patient meets both clinical and laboratory criteria.")
        final_answer = "Yes"
    else:
        print("Conclusion: The patient does not meet the full criteria for APS.")
        final_answer = "No"
        
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    diagnose_aps()