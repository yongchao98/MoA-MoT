def diagnose_antiphospholipid_syndrome():
    """
    Evaluates a patient's clinical and lab data against the 2006 Sydney revised
    classification criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    # Clinical Data
    vte_events = 3

    # Laboratory Data (T1 = 3 months ago, T2 = Today)
    labs_t1 = {
        'dRVVT_ratio': 1.44,
        'anticardiolipin_IgM': 32,
        'antiB2GP1_IgM': 41,
    }
    labs_t2 = {
        'dRVVT_ratio': 1.51,
        'anticardiolipin_IgM': 47,
        'antiB2GP1_IgM': 29,
    }
    # Normal Values for reference
    normals = {
        'dRVVT_ratio': 1.2,
        'aCL_IgM': 20,
        'aB2GP1_IgM': 20,
    }

    print("### Antiphospholipid Syndrome (APS) Diagnostic Evaluation ###")
    print("The diagnosis requires meeting at least one clinical AND one laboratory criterion.")
    print("-" * 60)

    # --- Step 1: Check Clinical Criteria ---
    print("1. Evaluating Clinical Criteria:")
    print(f"   Criterion: 1 or more episodes of venous, arterial, or small vessel thrombosis.")
    print(f"   Patient's History: The patient has experienced {vte_events} VTE events.")
    
    clinical_criterion_met = vte_events >= 1
    if clinical_criterion_met:
        print("   Result: Clinical criterion is MET.")
    else:
        print("   Result: Clinical criterion is NOT MET.")
    
    print("-" * 60)

    # --- Step 2: Check Laboratory Criteria ---
    print("2. Evaluating Laboratory Criteria:")
    print("   Criterion: A positive lab test on 2 or more occasions at least 12 weeks apart.")
    print("   The patient's tests are separated by 3 months, which meets the time requirement.")
    
    # Check for persistent Lupus Anticoagulant (LA) via dRVVT
    la_t1_positive = labs_t1['dRVVT_ratio'] > normals['dRVVT_ratio']
    la_t2_positive = labs_t2['dRVVT_ratio'] > normals['dRVVT_ratio']
    la_criterion_met = la_t1_positive and la_t2_positive
    print(f"\n   a) Lupus Anticoagulant (dRVVT):")
    print(f"      - 3 Months Ago: {labs_t1['dRVVT_ratio']} (Normal < {normals['dRVVT_ratio']}) -> Positive")
    print(f"      - Today: {labs_t2['dRVVT_ratio']} (Normal < {normals['dRVVT_ratio']}) -> Positive")
    print(f"      - Persistent? {'Yes' if la_criterion_met else 'No'}")

    # Check for persistent Anticardiolipin (aCL) IgM
    acl_t1_positive = labs_t1['anticardiolipin_IgM'] > normals['aCL_IgM']
    acl_t2_positive = labs_t2['anticardiolipin_IgM'] > normals['aCL_IgM']
    acl_criterion_met = acl_t1_positive and acl_t2_positive
    print(f"\n   b) Anticardiolipin IgM:")
    print(f"      - 3 Months Ago: {labs_t1['anticardiolipin_IgM']} UI/L (Normal < {normals['aCL_IgM']}) -> Positive")
    print(f"      - Today: {labs_t2['anticardiolipin_IgM']} UI/L (Normal < {normals['aCL_IgM']}) -> Positive")
    print(f"      - Persistent? {'Yes' if acl_criterion_met else 'No'}")

    # Check for persistent Anti-B2GP1 IgM
    b2gp1_t1_positive = labs_t1['antiB2GP1_IgM'] > normals['aB2GP1_IgM']
    b2gp1_t2_positive = labs_t2['antiB2GP1_IgM'] > normals['aB2GP1_IgM']
    b2gp1_criterion_met = b2gp1_t1_positive and b2gp1_t2_positive
    print(f"\n   c) Anti-B2GP1 IgM:")
    print(f"      - 3 Months Ago: {labs_t1['antiB2GP1_IgM']} UI/L (Normal < {normals['aB2GP1_IgM']}) -> Positive")
    print(f"      - Today: {labs_t2['antiB2GP1_IgM']} UI/L (Normal < {normals['aB2GP1_IgM']}) -> Positive")
    print(f"      - Persistent? {'Yes' if b2gp1_criterion_met else 'No'}")
    
    laboratory_criterion_met = la_criterion_met or acl_criterion_met or b2gp1_criterion_met
    print("\n   Result: At least one laboratory criterion is MET.")
    print("-" * 60)

    # --- Step 3: Final Conclusion ---
    print("3. Final Diagnosis:")
    if clinical_criterion_met and laboratory_criterion_met:
        final_answer = "Yes"
        print("   Conclusion: The patient meets BOTH clinical (thrombosis) and laboratory (persistent positive antibodies) criteria.")
    else:
        final_answer = "No"
        print("   Conclusion: The patient does NOT meet both required criteria for APS.")

    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    diagnose_antiphospholipid_syndrome()