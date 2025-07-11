def check_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the 2006 Sydney revised criteria.
    """

    # --- Patient Data ---
    # Clinical History
    vte_events = 3
    
    # Lab data from 3 months ago
    lab1 = {
        "anti_b2gp1_igm": 41, "anti_b2gp1_igg": 18,
        "anticardiolipin_igm": 32, "anticardiolipin_igg": 9,
        "ptt_la_ratio": 1.19, "drvv_ratio": 1.44
    }
    
    # Lab data from today (interval > 12 weeks)
    lab2 = {
        "anti_b2gp1_igm": 29, "anti_b2gp1_igg": 21,
        "anticardiolipin_igm": 47, "anticardiolipin_igg": 7,
        "ptt_la_ratio": 1.17, "drvv_ratio": 1.51
    }
    
    # Normal reference values
    normals = {
        "antibody_level": 20,
        "ptt_la_ratio": 1.18,
        "drvv_ratio": 1.2
    }

    print("Step 1: Evaluating Clinical Criteria for APS")
    print("-" * 40)
    
    clinical_criterion_met = False
    if vte_events >= 1:
        clinical_criterion_met = True
        print(f"Clinical Criterion: Met. The patient has a history of {vte_events} vascular thrombosis events.")
    else:
        print("Clinical Criterion: Not met.")
        
    print("\nStep 2: Evaluating Laboratory Criteria for APS")
    print("Criteria require at least one positive lab test on two occasions >12 weeks apart.")
    print("-" * 40)
    
    lab_criteria_met = False
    
    # Check 1: Lupus Anticoagulant (LA)
    la_positive_t1 = lab1["ptt_la_ratio"] > normals["ptt_la_ratio"] or lab1["drvv_ratio"] > normals["drvv_ratio"]
    la_positive_t2 = lab2["ptt_la_ratio"] > normals["ptt_la_ratio"] or lab2["drvv_ratio"] > normals["drvv_ratio"]
    
    print("Checking for persistent Lupus Anticoagulant (LA)...")
    print(f"  - Time 1: dRVVT ratio was {lab1['drvv_ratio']} (Normal < {normals['drvv_ratio']}). Positive: {lab1['drvv_ratio'] > normals['drvv_ratio']}.")
    print(f"  - Time 2: dRVVT ratio was {lab2['drvv_ratio']} (Normal < {normals['drvv_ratio']}). Positive: {lab2['drvv_ratio'] > normals['drvv_ratio']}.")
    if la_positive_t1 and la_positive_t2:
        lab_criteria_met = True
        print("  - Result: LA is persistently positive. Laboratory criterion MET.")
    else:
        print("  - Result: LA is not persistently positive.")

    # Check 2: Anticardiolipin (aCL) antibodies
    acl_positive_t1 = lab1["anticardiolipin_igm"] > normals["antibody_level"] or lab1["anticardiolipin_igg"] > normals["antibody_level"]
    acl_positive_t2 = lab2["anticardiolipin_igm"] > normals["antibody_level"] or lab2["anticardiolipin_igg"] > normals["antibody_level"]
    
    print("\nChecking for persistent Anticardiolipin (aCL) antibodies...")
    print(f"  - Time 1: aCL IgM was {lab1['anticardiolipin_igm']} UI/L (Normal < {normals['antibody_level']}). Positive: {lab1['anticardiolipin_igm'] > normals['antibody_level']}.")
    print(f"  - Time 2: aCL IgM was {lab2['anticardiolipin_igm']} UI/L (Normal < {normals['antibody_level']}). Positive: {lab2['anticardiolipin_igm'] > normals['antibody_level']}.")
    if acl_positive_t1 and acl_positive_t2:
        lab_criteria_met = True
        print("  - Result: aCL antibodies are persistently positive. Laboratory criterion MET.")
    else:
        print("  - Result: aCL antibodies are not persistently positive.")
        
    # Check 3: Anti-ß2GP1 antibodies
    b2gp1_positive_t1 = lab1["anti_b2gp1_igm"] > normals["antibody_level"] or lab1["anti_b2gp1_igg"] > normals["antibody_level"]
    b2gp1_positive_t2 = lab2["anti_b2gp1_igm"] > normals["antibody_level"] or lab2["anti_b2gp1_igg"] > normals["antibody_level"]

    print("\nChecking for persistent Anti-ß2GP1 antibodies...")
    print(f"  - Time 1: anti-ß2GP1 IgM was {lab1['anti_b2gp1_igm']} UI/L (Normal < {normals['antibody_level']}). Positive: {lab1['anti_b2gp1_igm'] > normals['antibody_level']}.")
    print(f"  - Time 2: anti-ß2GP1 IgM was {lab2['anti_b2gp1_igm']} UI/L and IgG was {lab2['anti_b2gp1_igg']} UI/L (Normal < {normals['antibody_level']}). Positive: {b2gp1_positive_t2}.")
    if b2gp1_positive_t1 and b2gp1_positive_t2:
        lab_criteria_met = True
        print("  - Result: Anti-ß2GP1 antibodies are persistently positive. Laboratory criterion MET.")
    else:
        print("  - Result: Anti-ß2GP1 antibodies are not persistently positive.")

    print("\nStep 3: Final Conclusion")
    print("-" * 40)

    if clinical_criterion_met and lab_criteria_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion AND at least one laboratory criterion.")
    else:
        final_answer = "No"
        print("The patient does not meet both the clinical and laboratory criteria for APS.")
        
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

check_aps_diagnosis()