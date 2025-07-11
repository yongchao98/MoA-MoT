import sys

def diagnose_antiphospholipid_syndrome():
    """
    Analyzes patient data to diagnose Antiphospholipid Syndrome (APS)
    based on the 2006 revised Sapporo (Sydney) classification criteria.
    """
    # --- Patient Data ---
    # Clinical History
    vte_events = 3

    # Lab results from 3 months ago (Test 1)
    # Using names T1_* for clarity
    T1_anticardiolipin_IgM = 32
    T1_antiB2GP1_IgM = 41
    T1_drvv_ratio = 1.44

    # Lab results from today (Test 2)
    # Using names T2_* for clarity
    T2_anticardiolipin_IgM = 47
    T2_antiB2GP1_IgM = 29
    T2_drvv_ratio = 1.51

    # --- Diagnostic Thresholds ---
    antibody_normal_titer = 20
    # Medium/high titer for aCL is often considered >40 units
    medium_titer_threshold = 40
    drvv_normal_ratio = 1.2
    
    print("Evaluating patient for Antiphospholipid Syndrome (APS)")
    print("Diagnosis requires at least 1 clinical AND 1 laboratory criterion.")
    print("="*50)

    # --- Step 1: Check Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria (Vascular Thrombosis)")
    
    # Check for one or more episodes of venous, arterial, or small-vessel thrombosis.
    clinical_criteria_met = vte_events >= 1
    
    print(f"Patient history shows {vte_events} VTE events.")
    print(f"Is the number of thrombotic events >= 1? {'Yes' if clinical_criteria_met else 'No'}")
    
    if clinical_criteria_met:
        print("Result: Clinical criterion for vascular thrombosis is MET.")
    else:
        print("Result: Clinical criterion is NOT MET.")
    print("="*50)

    # --- Step 2: Check Laboratory Criteria ---
    # Must be persistent (positive on 2+ occasions at least 12 weeks apart)
    print("Step 2: Evaluating Laboratory Criteria (Persistent Positivity)")
    print("The tests were performed >12 weeks apart, so persistence can be assessed.")

    # 2a. Check for persistent Lupus Anticoagulant (LA) via dRVVT
    print("\n--- Checking for persistent Lupus Anticoagulant (by dRVVT) ---")
    la_positive_T1 = T1_drvv_ratio > drvv_normal_ratio
    la_positive_T2 = T2_drvv_ratio > drvv_normal_ratio
    la_persistent = la_positive_T1 and la_positive_T2
    print(f"Test 1 dRVVT Ratio: {T1_drvv_ratio} (Normal < {drvv_normal_ratio}) -> Positive? {la_positive_T1}")
    print(f"Test 2 dRVVT Ratio: {T2_drvv_ratio} (Normal < {drvv_normal_ratio}) -> Positive? {la_positive_T2}")
    if la_persistent:
        print("Finding: Lupus Anticoagulant is PERSISTENTLY POSITIVE.")
    else:
        print("Finding: Lupus Anticoagulant is not persistently positive.")

    # 2b. Check for persistent Anticardiolipin (aCL) IgM antibodies
    print("\n--- Checking for persistent Anticardiolipin IgM (aCL IgM) ---")
    acl_positive_T1 = T1_anticardiolipin_IgM > antibody_normal_titer
    acl_positive_T2 = T2_anticardiolipin_IgM > antibody_normal_titer
    # Criterion requires medium-or-high titer, here defined as >40
    acl_medium_titer_present = T1_anticardiolipin_IgM > medium_titer_threshold or T2_anticardiolipin_IgM > medium_titer_threshold
    acl_persistent = acl_positive_T1 and acl_positive_T2 and acl_medium_titer_present
    print(f"Test 1 aCL IgM: {T1_anticardiolipin_IgM} (Normal < {antibody_normal_titer}) -> Positive? {acl_positive_T1}")
    print(f"Test 2 aCL IgM: {T2_anticardiolipin_IgM} (Normal < {antibody_normal_titer}) -> Positive? {acl_positive_T2}")
    if acl_persistent:
        print("Finding: aCL IgM is PERSISTENTLY POSITIVE at a medium/high titer.")
    else:
        print("Finding: aCL IgM is not persistently positive or not at medium/high titer.")

    # 2c. Check for persistent Anti-B2GP1 IgM antibodies
    print("\n--- Checking for persistent Anti-B2GP1 IgM ---")
    ab2gp1_positive_T1 = T1_antiB2GP1_IgM > antibody_normal_titer
    ab2gp1_positive_T2 = T2_antiB2GP1_IgM > antibody_normal_titer
    ab2gp1_persistent = ab2gp1_positive_T1 and ab2gp1_positive_T2
    print(f"Test 1 antiB2GP1 IgM: {T1_antiB2GP1_IgM} (Normal < {antibody_normal_titer}) -> Positive? {ab2gp1_positive_T1}")
    print(f"Test 2 antiB2GP1 IgM: {T2_antiB2GP1_IgM} (Normal < {antibody_normal_titer}) -> Positive? {ab2gp1_positive_T2}")
    if ab2gp1_persistent:
        print("Finding: Anti-B2GP1 IgM is PERSISTENTLY POSITIVE.")
    else:
        print("Finding: Anti-B2GP1 IgM is not persistently positive.")

    laboratory_criteria_met = la_persistent or acl_persistent or ab2gp1_persistent
    print("\n" + "-"*20)
    if laboratory_criteria_met:
        print("Result: At least one laboratory criterion is MET.")
    else:
        print("Result: Laboratory criteria are NOT MET.")
    print("="*50)
    
    # --- Step 3: Final Conclusion ---
    print("Step 3: Final Diagnosis")
    is_aps_diagnosed = clinical_criteria_met and laboratory_criteria_met
    
    print(f"Clinical Criteria Met? {'Yes' if clinical_criteria_met else 'No'}")
    print(f"Laboratory Criteria Met? {'Yes' if laboratory_criteria_met else 'No'}")
    
    if is_aps_diagnosed:
        final_answer = "Yes"
        print("\nCONCLUSION: The patient fulfills both clinical and laboratory criteria.")
        print("Diagnosis of Antiphospholipid Syndrome is confirmed.")
    else:
        final_answer = "No"
        print("\nCONCLUSION: The patient does not meet the full criteria for APS.")

    # Final Answer in requested format
    sys.stdout.write(f"\n<<<{final_answer}>>>\n")

if __name__ == '__main__':
    diagnose_antiphospholipid_syndrome()