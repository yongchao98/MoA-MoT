import sys

def check_aps_diagnosis():
    """
    Analyzes patient data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical History
    vte_events_count = 3 # Calf DVT, PE, another PE

    # Lab Test 1 (3 months ago)
    lab1 = {
        'dRVVT_ratio': 1.44,
        'dRVVT_norm': 1.2,
        'aCL_IgM': 32,
        'aCL_norm': 20,
        'aB2GP1_IgM': 41,
        'aB2GP1_norm': 20
    }

    # Lab Test 2 (Today - 12 weeks after test 1)
    lab2 = {
        'dRVVT_ratio': 1.51,
        'dRVVT_norm': 1.2,
        'aCL_IgM': 47,
        'aCL_norm': 20,
        'aB2GP1_IgM': 29,
        'aB2GP1_norm': 20
    }

    # Time between tests in weeks
    time_between_tests_weeks = 12

    # --- Step 1: Evaluate Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria (Sydney Criteria for APS)")
    print(f"A clinical criterion is at least one confirmed episode of venous or arterial thrombosis.")
    print(f"Patient's history includes {vte_events_count} VTE events.")
    
    clinical_criterion_met = vte_events_count >= 1
    if clinical_criterion_met:
        print("Result: Clinical criterion is MET.\n")
    else:
        print("Result: Clinical criterion is NOT MET.\n")

    # --- Step 2: Evaluate Laboratory Criteria ---
    print("Step 2: Evaluating Laboratory Criteria (Sydney Criteria for APS)")
    print("A laboratory criterion is the persistent presence of an antiphospholipid antibody on two occasions at least 12 weeks apart.\n")

    lab_criterion_met = False

    # Check for persistent Lupus Anticoagulant (using dRVVT)
    print("Checking for persistent Lupus Anticoagulant (LA)...")
    la_test1_positive = lab1['dRVVT_ratio'] > lab1['dRVVT_norm']
    la_test2_positive = lab2['dRVVT_ratio'] > lab2['dRVVT_norm']
    la_is_persistent = la_test1_positive and la_test2_positive and time_between_tests_weeks >= 12
    
    print(f"Test 1 dRVVT ratio: {lab1['dRVVT_ratio']} (Normal < {lab1['dRVVT_norm']}) -> Positive: {la_test1_positive}")
    print(f"Test 2 dRVVT ratio: {lab2['dRVVT_ratio']} (Normal < {lab2['dRVVT_norm']}) -> Positive: {la_test2_positive}")
    
    if la_is_persistent:
        lab_criterion_met = True
        print("Result: Persistent Lupus Anticoagulant DETECTED. Laboratory criterion is MET.\n")
    else:
        print("Result: No persistent Lupus Anticoagulant detected.\n")

    # Check for persistent anticardiolipin antibodies (IgM)
    if not lab_criterion_met:
        print("Checking for persistent Anticardiolipin (aCL) IgM...")
        acl_test1_positive = lab1['aCL_IgM'] > lab1['aCL_norm']
        acl_test2_positive = lab2['aCL_IgM'] > lab2['aCL_norm']
        acl_is_persistent = acl_test1_positive and acl_test2_positive and time_between_tests_weeks >= 12

        print(f"Test 1 aCL IgM: {lab1['aCL_IgM']} UI/L (Normal < {lab1['aCL_norm']}) -> Positive: {acl_test1_positive}")
        print(f"Test 2 aCL IgM: {lab2['aCL_IgM']} UI/L (Normal < {lab2['aCL_norm']}) -> Positive: {acl_test2_positive}")

        if acl_is_persistent:
            lab_criterion_met = True
            print("Result: Persistent aCL IgM DETECTED. Laboratory criterion is MET.\n")
        else:
            print("Result: No persistent aCL IgM detected.\n")

    # Check for persistent anti-B2GP1 antibodies (IgM)
    if not lab_criterion_met:
        print("Checking for persistent Anti-ß2GP1 IgM...")
        ab2gp1_test1_positive = lab1['aB2GP1_IgM'] > lab1['aB2GP1_norm']
        ab2gp1_test2_positive = lab2['aB2GP1_IgM'] > lab2['aB2GP1_norm']
        ab2gp1_is_persistent = ab2gp1_test1_positive and ab2gp1_test2_positive and time_between_tests_weeks >= 12

        print(f"Test 1 aB2GP1 IgM: {lab1['aB2GP1_IgM']} UI/L (Normal < {lab1['aB2GP1_norm']}) -> Positive: {ab2gp1_test1_positive}")
        print(f"Test 2 aB2GP1 IgM: {lab2['aB2GP1_IgM']} UI/L (Normal < {lab2['aB2GP1_norm']}) -> Positive: {ab2gp1_test2_positive}")

        if ab2gp1_is_persistent:
            lab_criterion_met = True
            print("Result: Persistent Anti-ß2GP1 IgM DETECTED. Laboratory criterion is MET.\n")
        else:
            print("Result: No persistent Anti-ß2GP1 IgM detected.\n")

    # --- Step 3: Final Conclusion ---
    print("Step 3: Final Conclusion")
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {lab_criterion_met}")
    
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
    else:
        final_answer = "No"

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    
    # Writing the final answer to stdout for grading.
    # The format <<<ANSWER>>> is a specified requirement.
    sys.stdout.write(f'<<<{final_answer}>>>')

check_aps_diagnosis()