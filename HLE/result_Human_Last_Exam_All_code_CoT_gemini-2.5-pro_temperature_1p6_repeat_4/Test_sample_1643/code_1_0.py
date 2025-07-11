def check_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the 2006 revised Sapporo criteria.
    """
    # --- Patient Data ---
    # Clinical History
    vte_events = 3

    # Lab tests 3 months ago
    lab_1 = {
        'dRVVT_ratio': 1.44,
        'aCL_IgM': 32,
        'aCL_IgG': 9,
        'aB2GP1_IgM': 41,
        'aB2GP1_IgG': 18
    }

    # Lab tests today (interval > 12 weeks)
    lab_2 = {
        'dRVVT_ratio': 1.51,
        'aCL_IgM': 47,
        'aCL_IgG': 7,
        'aB2GP1_IgM': 29,
        'aB2GP1_IgG': 21
    }

    # --- Criteria Thresholds ---
    vte_threshold = 1
    dRVVT_threshold = 1.2
    antibody_threshold = 20 # UI/L, normal is < 20

    print("Step 1: Evaluating Clinical Criteria for APS")
    print("------------------------------------------")
    # Check for Vascular Thrombosis
    clinical_criterion_met = vte_events >= vte_threshold
    print(f"Patient has {vte_events} documented VTE events.")
    print(f"The criterion requires at least {vte_threshold} event.")
    print(f"Is {vte_events} >= {vte_threshold}? {'Yes' if clinical_criterion_met else 'No'}.")
    print(f"Conclusion: Clinical criterion is {'MET' if clinical_criterion_met else 'NOT MET'}.\n")

    print("Step 2: Evaluating Laboratory Criteria for APS")
    print("---------------------------------------------")
    print("Note: Criteria require persistence on two tests at least 12 weeks apart.\n")

    # Check for persistent Lupus Anticoagulant (LA) via dRVVT
    la_positive_1 = lab_1['dRVVT_ratio'] > dRVVT_threshold
    la_positive_2 = lab_2['dRVVT_ratio'] > dRVVT_threshold
    la_criterion_met = la_positive_1 and la_positive_2
    print("Criterion A: Persistent Lupus Anticoagulant")
    print(f" - Test 1 (3 months ago): dRVVT ratio is {lab_1['dRVVT_ratio']}. Is this > {dRVVT_threshold}? {'Yes' if la_positive_1 else 'No'}.")
    print(f" - Test 2 (Today): dRVVT ratio is {lab_2['dRVVT_ratio']}. Is this > {dRVVT_threshold}? {'Yes' if la_positive_2 else 'No'}.")
    print(f"Conclusion: Lupus Anticoagulant criterion is {'MET' if la_criterion_met else 'NOT MET'}.\n")

    # Check for persistent Anticardiolipin (aCL) IgM/IgG antibodies
    acl_igm_positive_1 = lab_1['aCL_IgM'] > antibody_threshold
    acl_igm_positive_2 = lab_2['aCL_IgM'] > antibody_threshold
    persistent_acl_igm = acl_igm_positive_1 and acl_igm_positive_2
    
    acl_igg_positive_1 = lab_1['aCL_IgG'] > antibody_threshold
    acl_igg_positive_2 = lab_2['aCL_IgG'] > antibody_threshold
    persistent_acl_igg = acl_igg_positive_1 and acl_igg_positive_2
    
    acl_criterion_met = persistent_acl_igm or persistent_acl_igg
    print("Criterion B: Persistent Anticardiolipin Antibodies (IgM or IgG)")
    print(f" - Test 1 (3 months ago): aCL IgM is {lab_1['aCL_IgM']} (> {antibody_threshold}? {'Yes' if acl_igm_positive_1 else 'No'}); aCL IgG is {lab_1['aCL_IgG']} (> {antibody_threshold}? {'Yes' if acl_igg_positive_1 else 'No'}).")
    print(f" - Test 2 (Today): aCL IgM is {lab_2['aCL_IgM']} (> {antibody_threshold}? {'Yes' if acl_igm_positive_2 else 'No'}); aCL IgG is {lab_2['aCL_IgG']} (> {antibody_threshold}? {'Yes' if acl_igg_positive_2 else 'No'}).")
    print(f"Conclusion: Anticardiolipin Antibody criterion is {'MET' if acl_criterion_met else 'NOT MET'}.\n")

    # Check for persistent Anti-B2GP1 IgM/IgG antibodies
    ab2gp1_igm_positive_1 = lab_1['aB2GP1_IgM'] > antibody_threshold
    ab2gp1_igm_positive_2 = lab_2['aB2GP1_IgM'] > antibody_threshold
    persistent_ab2gp1_igm = ab2gp1_igm_positive_1 and ab2gp1_igm_positive_2

    ab2gp1_igg_positive_1 = lab_1['aB2GP1_IgG'] > antibody_threshold
    ab2gp1_igg_positive_2 = lab_2['aB2GP1_IgG'] > antibody_threshold
    persistent_ab2gp1_igg = ab2gp1_igg_positive_1 and ab2gp1_igg_positive_2

    ab2gp1_criterion_met = persistent_ab2gp1_igm or persistent_ab2gp1_igg
    print("Criterion C: Persistent Anti-β2 Glycoprotein-I Antibodies (IgM or IgG)")
    print(f" - Test 1 (3 months ago): aβ2GPI IgM is {lab_1['aB2GP1_IgM']} (> {antibody_threshold}? {'Yes' if ab2gp1_igm_positive_1 else 'No'}); aβ2GPI IgG is {lab_1['aB2GP1_IgG']} (> {antibody_threshold}? {'Yes' if ab2gp1_igg_positive_1 else 'No'}).")
    print(f" - Test 2 (Today): aβ2GPI IgM is {lab_2['aB2GP1_IgM']} (> {antibody_threshold}? {'Yes' if ab2gp1_igm_positive_2 else 'No'}); aβ2GPI IgG is {lab_2['aB2GP1_IgG']} (> {antibody_threshold}? {'Yes' if ab2gp1_igg_positive_2 else 'No'}).")
    print(f"Conclusion: Anti-β2 Glycoprotein-I Antibody criterion is {'MET' if ab2gp1_criterion_met else 'NOT MET'}.\n")

    # Overall laboratory criteria check
    lab_criteria_met = la_criterion_met or acl_criterion_met or ab2gp1_criterion_met
    print("Final check of Laboratory Criteria:")
    print(f"At least one laboratory criterion must be met. Is this condition satisfied? {'Yes' if lab_criteria_met else 'No'}.")
    print(f"Overall Laboratory Criteria: {'MET' if lab_criteria_met else 'NOT MET'}.\n")


    print("Step 3: Final Diagnosis")
    print("------------------------")
    final_diagnosis = clinical_criterion_met and lab_criteria_met
    print(f"To diagnose APS, we need: (Clinical Criteria MET) AND (Laboratory Criteria MET).")
    print(f"In this case: ({'MET' if clinical_criterion_met else 'NOT MET'}) AND ({'MET' if lab_criteria_met else 'NOT MET'}).")
    
    if final_diagnosis:
        final_answer = "Yes"
    else:
        final_answer = "No"

    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    check_aps_diagnosis()