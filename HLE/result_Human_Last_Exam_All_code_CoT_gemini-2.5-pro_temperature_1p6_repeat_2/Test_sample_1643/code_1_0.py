import sys

# Redirect print to a string buffer to capture output for the final response format if needed
# This is a good practice for complex outputs, but for this task, direct printing is fine.

def diagnose_aps():
    """
    Evaluates patient data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical History
    vte_events_count = 3

    # Laboratory Data - Test 1 (3 months ago)
    lab1 = {
        'dRVVT_ratio': 1.44, 'dRVVT_norm': 1.2,
        'aCL_IgM': 32, 'aCL_norm': 20,
        'aCL_IgG': 9,
        'ab2GPI_IgM': 41, 'ab2GPI_norm': 20,
        'ab2GPI_IgG': 18
    }

    # Laboratory Data - Test 2 (Today)
    lab2 = {
        'dRVVT_ratio': 1.51, 'dRVVT_norm': 1.2,
        'aCL_IgM': 47, 'aCL_norm': 20,
        'aCL_IgG': 7,
        'ab2GPI_IgM': 29, 'ab2GPI_norm': 20,
        'ab2GPI_IgG': 21
    }

    # Time between tests in months
    time_between_tests_months = 3
    
    # --- Evaluation Logic ---
    print("Evaluating patient for Antiphospholipid Syndrome (APS) using the Sydney criteria.")
    print("-" * 60)

    # 1. Check Clinical Criteria
    print("Step 1: Analyzing Clinical Criteria...")
    print(f"The patient has a history of {vte_events_count} VTE events.")
    clinical_criterion_met = vte_events_count >= 1
    if clinical_criterion_met:
        print("Result: Patient meets the clinical criterion for vascular thrombosis (>=1 event).")
    else:
        print("Result: Patient does not meet the clinical criterion for vascular thrombosis.")
    print("-" * 60)

    # 2. Check Laboratory Criteria
    print("Step 2: Analyzing Laboratory Criteria (must be persistent over >= 12 weeks)...")
    print(f"Lab tests were performed {time_between_tests_months} months apart, meeting the time requirement.")
    
    lab_criterion_met = False
    
    # Check for persistent Lupus Anticoagulant (LA) - using dRVVT as the primary indicator
    la_test1_positive = lab1['dRVVT_ratio'] > lab1['dRVVT_norm']
    la_test2_positive = lab2['dRVVT_ratio'] > lab2['dRVVT_norm']
    persistent_la = la_test1_positive and la_test2_positive
    print(f"Lupus Anticoagulant (dRVVT):")
    print(f"  - Test 1: {lab1['dRVVT_ratio']} (Normal < {lab1['dRVVT_norm']}) -> {'Positive' if la_test1_positive else 'Negative'}")
    print(f"  - Test 2: {lab2['dRVVT_ratio']} (Normal < {lab2['dRVVT_norm']}) -> {'Positive' if la_test2_positive else 'Negative'}")
    if persistent_la:
        print("  -> Persistent LA found. Laboratory criterion met.")
        lab_criterion_met = True
    else:
        print("  -> LA is not persistently positive.")
        
    # Check for persistent anticardiolipin (aCL) antibodies
    acl_test1_positive = lab1['aCL_IgM'] > lab1['aCL_norm'] or lab1['aCL_IgG'] > lab1['aCL_norm']
    acl_test2_positive = lab2['aCL_IgM'] > lab2['aCL_norm'] or lab2['aCL_IgG'] > lab2['aCL_norm']
    persistent_acl = acl_test1_positive and acl_test2_positive
    print(f"Anticardiolipin Antibodies (aCL):")
    print(f"  - Test 1 (IgM): {lab1['aCL_IgM']} (Normal < {lab1['aCL_norm']}) -> {'Positive' if lab1['aCL_IgM'] > lab1['aCL_norm'] else 'Negative'}")
    print(f"  - Test 2 (IgM): {lab2['aCL_IgM']} (Normal < {lab2['aCL_norm']}) -> {'Positive' if lab2['aCL_IgM'] > lab2['aCL_norm'] else 'Negative'}")
    if persistent_acl:
        print("  -> Persistent aCL antibodies found. Laboratory criterion met.")
        if not lab_criterion_met: lab_criterion_met = True # Set if not already set by LA
    else:
        print("  -> aCL antibodies are not persistently positive.")
        
    # Check for persistent anti-beta2-glycoprotein-I (aB2GPI) antibodies
    ab2gpi_test1_positive = lab1['ab2GPI_IgM'] > lab1['ab2GPI_norm'] or lab1['ab2GPI_IgG'] > lab1['ab2GPI_norm']
    ab2gpi_test2_positive = lab2['ab2GPI_IgM'] > lab2['ab2GPI_norm'] or lab2['ab2GPI_IgG'] > lab2['ab2GPI_norm']
    persistent_ab2gpi = ab2gpi_test1_positive and ab2gpi_test2_positive
    print(f"Anti-B2GP1 Antibodies (aB2GPI):")
    print(f"  - Test 1 (IgM): {lab1['ab2GPI_IgM']} (Normal < {lab1['ab2GPI_norm']}) -> {'Positive' if lab1['ab2GPI_IgM'] > lab1['ab2GPI_norm'] else 'Negative'}")
    print(f"  - Test 2 (IgM/IgG): IgM {lab2['ab2GPI_IgM']}, IgG {lab2['ab2GPI_IgG']} (Normal < {lab2['ab2GPI_norm']}) -> {'Positive' if ab2gpi_test2_positive else 'Negative'}")
    if persistent_ab2gpi:
        print("  -> Persistent aB2GPI antibodies found. Laboratory criterion met.")
        if not lab_criterion_met: lab_criterion_met = True
    else:
        print("  -> aB2GPI antibodies are not persistently positive.")

    print("\nResult: Patient meets at least one laboratory criterion.") if lab_criterion_met else print("\nResult: Patient does not meet any laboratory criterion.")
    print("-" * 60)
    
    # 3. Final Conclusion
    print("Step 3: Final Conclusion...")
    has_aps = clinical_criterion_met and lab_criterion_met
    
    if has_aps:
        print("The patient meets both clinical and laboratory criteria for Antiphospholipid Syndrome.")
        final_answer = "Yes"
    else:
        print("The patient does not meet the full criteria for Antiphospholipid Syndrome.")
        final_answer = "No"
        
    return final_answer

if __name__ == '__main__':
    answer = diagnose_aps()
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?\n<<<{answer}>>>")