def diagnose_aps():
    """
    Analyzes patient data to determine if they meet the criteria for Antiphospholipid Syndrome (APS).
    """
    # --- Patient Data ---
    # Clinical
    vte_events = 3
    
    # Lab Test 1 (3 months ago)
    lab1 = {
        "dRVVT_ratio": 1.44,
        "aCL_IgM": 32,
        "aCL_IgG": 9,
        "aB2GP1_IgM": 41,
        "aB2GP1_IgG": 18
    }
    
    # Lab Test 2 (Today, >12 weeks later)
    lab2 = {
        "dRVVT_ratio": 1.51,
        "aCL_IgM": 47,
        "aCL_IgG": 7,
        "aB2GP1_IgM": 29,
        "aB2GP1_IgG": 21
    }
    
    # Normal Values
    normals = {
        "dRVVT_ratio_max": 1.2,
        "antibody_max": 20
    }

    print("Analyzing patient data for Antiphospholipid Syndrome (APS) diagnosis...")
    print("-" * 30)

    # --- Step 1: Check Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria")
    clinical_criteria_met = False
    if vte_events >= 1:
        clinical_criteria_met = True
        print(f"-> Finding: Patient has a history of {vte_events} venous thromboembolism (VTE) events.")
        print("-> Conclusion: The clinical criterion (vascular thrombosis) is MET.")
    else:
        print("-> Conclusion: The clinical criterion is NOT MET.")
    
    print("-" * 30)

    # --- Step 2: Check Laboratory Criteria ---
    print("Step 2: Evaluating Laboratory Criteria (Persistence over >12 weeks)")
    lab_criteria_met = False
    
    # Check for persistent Lupus Anticoagulant (using dRVVT)
    la_positive_test1 = lab1["dRVVT_ratio"] > normals["dRVVT_ratio_max"]
    la_positive_test2 = lab2["dRVVT_ratio"] > normals["dRVVT_ratio_max"]
    persistent_la = la_positive_test1 and la_positive_test2
    
    print(f"-> Checking for persistent Lupus Anticoagulant (dRVVT > {normals['dRVVT_ratio_max']}):")
    print(f"  - Test 1 (3 months ago): dRVVT = {lab1['dRVVT_ratio']} (Positive: {la_positive_test1})")
    print(f"  - Test 2 (Today): dRVVT = {lab2['dRVVT_ratio']} (Positive: {la_positive_test2})")
    if persistent_la:
        lab_criteria_met = True
        print("  - Result: Persistent Lupus Anticoagulant DETECTED.")
    else:
        print("  - Result: Persistent Lupus Anticoagulant NOT detected.")

    # Check for persistent anticardiolipin (aCL) antibodies
    acl_positive_test1 = lab1["aCL_IgM"] > normals["antibody_max"] or lab1["aCL_IgG"] > normals["antibody_max"]
    acl_positive_test2 = lab2["aCL_IgM"] > normals["antibody_max"] or lab2["aCL_IgG"] > normals["antibody_max"]
    persistent_acl = acl_positive_test1 and acl_positive_test2

    print(f"-> Checking for persistent anticardiolipin antibodies (> {normals['antibody_max']} UI/L):")
    print(f"  - Test 1 (3 months ago): aCL IgM = {lab1['aCL_IgM']}, aCL IgG = {lab1['aCL_IgG']} (Positive: {acl_positive_test1})")
    print(f"  - Test 2 (Today): aCL IgM = {lab2['aCL_IgM']}, aCL IgG = {lab2['aCL_IgG']} (Positive: {acl_positive_test2})")
    if persistent_acl:
        lab_criteria_met = True
        print("  - Result: Persistent anticardiolipin antibodies DETECTED.")
    else:
        print("  - Result: Persistent anticardiolipin antibodies NOT detected.")

    # Check for persistent anti-B2GP1 antibodies
    ab2gp1_positive_test1 = lab1["aB2GP1_IgM"] > normals["antibody_max"] or lab1["aB2GP1_IgG"] > normals["antibody_max"]
    ab2gp1_positive_test2 = lab2["aB2GP1_IgM"] > normals["antibody_max"] or lab2["aB2GP1_IgG"] > normals["antibody_max"]
    persistent_ab2gp1 = ab2gp1_positive_test1 and ab2gp1_positive_test2

    print(f"-> Checking for persistent anti-B2GP1 antibodies (> {normals['antibody_max']} UI/L):")
    print(f"  - Test 1 (3 months ago): aB2GP1 IgM = {lab1['aB2GP1_IgM']}, aB2GP1 IgG = {lab1['aB2GP1_IgG']} (Positive: {ab2gp1_positive_test1})")
    print(f"  - Test 2 (Today): aB2GP1 IgM = {lab2['aB2GP1_IgM']}, aB2GP1 IgG = {lab2['aB2GP1_IgG']} (Positive: {ab2gp1_positive_test2})")
    if persistent_ab2gp1:
        lab_criteria_met = True
        print("  - Result: Persistent anti-B2GP1 antibodies DETECTED.")
    else:
        print("  - Result: Persistent anti-B2GP1 antibodies NOT detected.")

    if lab_criteria_met:
        print("-> Conclusion: At least one laboratory criterion is MET.")
    else:
        print("-> Conclusion: Laboratory criteria are NOT MET.")
        
    print("-" * 30)

    # --- Step 3: Final Diagnosis ---
    print("Step 3: Final Diagnosis")
    print(f"Clinical Criteria Met: {clinical_criteria_met}")
    print(f"Laboratory Criteria Met: {lab_criteria_met}")
    
    if clinical_criteria_met and lab_criteria_met:
        final_answer = "Yes"
        print("\nConclusion: The patient meets both clinical and laboratory criteria for Antiphospholipid Syndrome.")
    else:
        final_answer = "No"
        print("\nConclusion: The patient does not meet the full criteria for Antiphospholipid Syndrome.")
        
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    diagnose_aps()