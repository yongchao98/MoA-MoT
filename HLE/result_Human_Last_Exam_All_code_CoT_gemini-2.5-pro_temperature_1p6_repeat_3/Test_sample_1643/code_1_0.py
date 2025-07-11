def diagnose_aps():
    """
    Analyzes patient data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical History
    vte_events = 3
    
    # Lab tests from 3 months ago (Time 1)
    lab_t1 = {
        'ab2gp1_igm': 41,  # N < 20 UI/L
        'ab2gp1_igg': 18,  # N < 20 UI/L
        'ac_igm': 32,      # N < 20 UI/L
        'ac_igg': 9,       # N < 20 UI/L
        'drvvt_ratio': 1.44 # N < 1.2
    }

    # Lab tests from today (Time 2)
    lab_t2 = {
        'ab2gp1_igm': 29,  # N < 20 UI/L
        'ab2gp1_igg': 21,  # N < 20 UI/L
        'ac_igm': 47,      # N < 20 UI/L
        'ac_igg': 7,       # N < 20 UI/L
        'drvvt_ratio': 1.51 # N < 1.2
    }
    
    # --- Analysis ---
    print("Step 1: Evaluating Clinical Criteria for APS")
    
    clinical_criterion_met = False
    if vte_events >= 1:
        clinical_criterion_met = True
        print(f"Patient has a history of {vte_events} VTE events.")
        print("-> Clinical Criterion for Vascular Thrombosis: MET\n")
    else:
        print("-> Clinical Criterion for Vascular Thrombosis: NOT MET\n")

    print("Step 2: Evaluating Laboratory Criteria for APS (Persistence over >12 weeks)")
    print("The lab tests were taken 3 months apart, which satisfies the requirement.\n")

    lab_criterion_met = False
    
    # Check for persistent Lupus Anticoagulant (LA) via dRVVT
    # Note: dRVVT can be falsely positive due to Rivaroxaban, but we will evaluate per the numbers given.
    print("Checking for persistent Lupus Anticoagulant (dRVVT > 1.2):")
    la_t1_positive = lab_t1['drvvt_ratio'] > 1.2
    la_t2_positive = lab_t2['drvvt_ratio'] > 1.2
    print(f"Time 1: dRVVT Ratio = {lab_t1['drvvt_ratio']}. Positive: {la_t1_positive}")
    print(f"Time 2: dRVVT Ratio = {lab_t2['drvvt_ratio']}. Positive: {la_t2_positive}")
    if la_t1_positive and la_t2_positive:
        lab_criterion_met = True
        print("-> Laboratory Criterion for Lupus Anticoagulant: MET\n")
    else:
        print("-> Laboratory Criterion for Lupus Anticoagulant: NOT MET\n")

    # Check for persistent Anticardiolipin (aCL) antibodies
    print("Checking for persistent Anticardiolipin (aCL) antibodies (IgM or IgG):")
    # IgM check
    acl_igm_persistent = lab_t1['ac_igm'] > 20 and lab_t2['ac_igm'] > 20
    acl_igm_medium_high_titer = lab_t1['ac_igm'] > 40 or lab_t2['ac_igm'] > 40
    print(f"aCL IgM at Time 1 = {lab_t1['ac_igm']}; at Time 2 = {lab_t2['ac_igm']}. Both are > 20 (Normal < 20).")
    if acl_igm_persistent and acl_igm_medium_high_titer:
        lab_criterion_met = True
        print(f"At least one titer is > 40 (Time 2: {lab_t2['ac_igm']}).")
        print("-> Laboratory Criterion for aCL IgM: MET\n")
    else:
        print("-> Laboratory Criterion for aCL IgM: NOT MET\n")
        
    # Check for persistent Anti-ß2GP1 antibodies
    print("Checking for persistent Anti-ß2GP1 antibodies (IgM or IgG):")
    # IgM check
    ab2gp1_igm_persistent = lab_t1['ab2gp1_igm'] > 20 and lab_t2['ab2gp1_igm'] > 20
    ab2gp1_igm_medium_high_titer = lab_t1['ab2gp1_igm'] > 40 or lab_t2['ab2gp1_igm'] > 40
    print(f"Anti-ß2GP1 IgM at Time 1 = {lab_t1['ab2gp1_igm']}; at Time 2 = {lab_t2['ab2gp1_igm']}. Both are > 20 (Normal < 20).")
    if ab2gp1_igm_persistent and ab2gp1_igm_medium_high_titer:
        lab_criterion_met = True
        print(f"At least one titer is > 40 (Time 1: {lab_t1['ab2gp1_igm']}).")
        print("-> Laboratory Criterion for Anti-ß2GP1 IgM: MET\n")
    else:
        print("-> Laboratory Criterion for Anti-ß2GP1 IgM: NOT MET\n")


    # --- Final Conclusion ---
    print("Step 3: Final Conclusion")
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {lab_criterion_met}")

    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("\nThe patient meets at least one clinical and one laboratory criterion for Antiphospholipid Syndrome.")
    else:
        final_answer = "No"
        print("\nThe patient does not meet the criteria for Antiphospholipid Syndrome.")
        
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    diagnose_aps()