import sys

def analyze_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) and prints the reasoning.
    """
    # --- Patient Data ---
    # Clinical
    vte_events = 3
    
    # Lab Test 1 (3 months ago)
    lab_test_1 = {
        "anti_b2gp1_igm": 41,
        "anticardiolipin_igm": 32,
        "dRVVT_ratio": 1.44
    }
    
    # Lab Test 2 (Today)
    lab_test_2 = {
        "anti_b2gp1_igm": 29,
        "anti_b2gp1_igg": 21,
        "anticardiolipin_igm": 47,
        "dRVVT_ratio": 1.51
    }
    
    # Normal reference values
    normal_ranges = {
        "antibody_titer": 20, # For aCL and anti-b2gp1
        "dRVVT_ratio": 1.2
    }
    
    # Initialize criteria flags
    clinical_criterion_met = False
    lab_criterion_met = False
    
    print("Evaluating patient for Antiphospholipid Syndrome (APS)...")
    print("-" * 50)
    
    # 1. Evaluate Clinical Criteria
    print("1. Clinical Criteria Evaluation (Requires at least 1 VTE event):")
    print(f"   - Patient history shows {vte_events} VTE events.")
    if vte_events >= 1:
        clinical_criterion_met = True
        print("   -> FINDING: The clinical criterion for vascular thrombosis is MET.")
    else:
        print("   -> FINDING: The clinical criterion is NOT met.")
        
    print("-" * 50)
    
    # 2. Evaluate Laboratory Criteria (Requires persistent positive test > 12 weeks)
    print("2. Laboratory Criteria Evaluation (Tests are 3 months apart, satisfying >12 week interval):")
    
    # a) Lupus Anticoagulant (LA) via dRVVT
    dRVVT_1 = lab_test_1['dRVVT_ratio']
    dRVVT_2 = lab_test_2['dRVVT_ratio']
    dRVVT_norm = normal_ranges['dRVVT_ratio']
    la_positive = dRVVT_1 > dRVVT_norm and dRVVT_2 > dRVVT_norm
    print(f"   a) Lupus Anticoagulant (LA) positivity (dRVVT > {dRVVT_norm}):")
    print(f"      - Test 1 dRVVT Ratio: {dRVVT_1}")
    print(f"      - Test 2 dRVVT Ratio: {dRVVT_2}")
    if la_positive:
        print("      -> FINDING: LA is persistently positive.")
        
    # b) Anticardiolipin (aCL) IgM
    acl_1 = lab_test_1['anticardiolipin_igm']
    acl_2 = lab_test_2['anticardiolipin_igm']
    acl_norm = normal_ranges['antibody_titer']
    acl_positive = acl_1 > acl_norm and acl_2 > acl_norm
    print(f"\n   b) Anticardiolipin IgM positivity (aCL > {acl_norm} UI/L):")
    print(f"      - Test 1 aCL IgM: {acl_1} UI/L")
    print(f"      - Test 2 aCL IgM: {acl_2} UI/L")
    if acl_positive:
        print("      -> FINDING: Anticardiolipin IgM is persistently positive.")
        
    # c) Anti-β2GP1
    b2gp1_igm1 = lab_test_1['anti_b2gp1_igm']
    b2gp1_igm2 = lab_test_2['anti_b2gp1_igm']
    b2gp1_igg2 = lab_test_2['anti_b2gp1_igg']
    b2gp1_norm = normal_ranges['antibody_titer']
    b2gp1_positive_t1 = b2gp1_igm1 > b2gp1_norm
    b2gp1_positive_t2 = b2gp1_igm2 > b2gp1_norm or b2gp1_igg2 > b2gp1_norm
    b2gp1_positive = b2gp1_positive_t1 and b2gp1_positive_t2
    print(f"\n   c) Anti-β2GP1 positivity (anti-β2GP1 > {b2gp1_norm} UI/L):")
    print(f"      - Test 1 anti-β2GP1 IgM: {b2gp1_igm1} UI/L")
    print(f"      - Test 2 anti-β2GP1 IgM: {b2gp1_igm2} UI/L; IgG: {b2gp1_igg2} UI/L")
    if b2gp1_positive:
        print("      -> FINDING: Anti-β2GP1 antibody is persistently positive.")
        
    if la_positive or acl_positive or b2gp1_positive:
        lab_criterion_met = True
        print("\n   -> FINDING: At least one laboratory criterion is MET.")
    else:
        print("\n   -> FINDING: The laboratory criterion is NOT met.")
        
    print("-" * 50)
    
    # 3. Final Conclusion
    print("3. Final Diagnosis:")
    final_answer = "No"
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets both clinical (recurrent VTE) and laboratory (persistent positive antiphospholipid antibodies) criteria.")
        print("\nCONCLUSION: The patient categorizes as having antiphospholipid syndrome.")
    else:
        print("The patient does not meet both the required clinical and laboratory criteria.")
        print("\nCONCLUSION: The patient does not categorize as having antiphospholipid syndrome.")
        
    # This is a special instruction for the final answer format
    sys.stdout.write(f"\n<<<{final_answer}>>>\n")

# Run the analysis
analyze_aps_diagnosis()