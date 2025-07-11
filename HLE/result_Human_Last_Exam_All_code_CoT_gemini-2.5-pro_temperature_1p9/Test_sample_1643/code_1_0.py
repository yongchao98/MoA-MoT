def solve_aps_diagnosis():
    """
    Analyzes patient data against the 2006 Sydney revised criteria for Antiphospholipid Syndrome (APS)
    and prints a step-by-step diagnostic conclusion.
    """
    # --- Patient Data ---
    patient_data = {
        'vte_events': 3,
        'medications': ['DESOGESTREL 75µg', 'RIVAROXABAN 20mg'],
        'labs': [
            {
                'time': '3 months ago',
                'aCL_IgM': 32, 'aCL_IgG': 9,
                'ab2gp1_IgM': 41, 'ab2gp1_IgG': 18,
                'drvvtratio': 1.44
            },
            {
                'time': 'Today',
                'aCL_IgM': 47, 'aCL_IgG': 7,
                'ab2gp1_IgM': 29, 'ab2gp1_IgG': 21,
                'drvvtratio': 1.51
            }
        ]
    }

    # --- APS Diagnostic Criteria Thresholds (Sydney 2006) ---
    criteria = {
        'vte_min': 1,
        'norm_acl_ab2gp1': 20,
        'titer_acl_ab2gp1': 40, # Medium/high titer cutoff (or >99th percentile)
        'norm_drvvtratio': 1.2
    }

    print("Analyzing patient case for Antiphospholipid Syndrome (APS) using the 2006 Sydney revised criteria.")
    print("-" * 30)

    # --- 1. Check Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria")
    clinical_met = False
    vte_count = patient_data['vte_events']
    print(f"The patient has a history of {vte_count} venous thromboembolism (VTE) events.")
    if vte_count >= criteria['vte_min']:
        clinical_met = True
        print("-> Clinical criterion for Vascular Thrombosis is MET.")
    else:
        print("-> Clinical criterion for Vascular Thrombosis is NOT MET.")
    print("-" * 30)

    # --- 2. Check Laboratory Criteria ---
    print("Step 2: Evaluating Laboratory Criteria")
    print("Note: Criteria require persistence on two occasions at least 12 weeks apart.")
    
    lab_met = False
    la_met = False
    acl_met = False
    ab2gp1_met = False
    
    labs1 = patient_data['labs'][0]
    labs2 = patient_data['labs'][1]

    # 2a. Lupus Anticoagulant (LA)
    print("\n* Checking for Lupus Anticoagulant (LA):")
    la_test_1_pos = labs1['drvvtratio'] > criteria['norm_drvvtratio']
    la_test_2_pos = labs2['drvvtratio'] > criteria['norm_drvvtratio']
    
    print(f"  - 3 months ago: dRVVT ratio was {labs1['drvvtratio']} (Normal < {criteria['norm_drvvtratio']}). Result: {'Positive' if la_test_1_pos else 'Negative'}.")
    print(f"  - Today: dRVVT ratio is {labs2['drvvtratio']} (Normal < {criteria['norm_drvvtratio']}). Result: {'Positive' if la_test_2_pos else 'Negative'}.")

    if la_test_1_pos and la_test_2_pos:
        la_met = True
        print("  -> Finding: Persistent positive LA detected based on dRVVT.")
        print("     -> LA laboratory criterion is MET.")
        if 'RIVAROXABAN 20mg' in patient_data['medications']:
             print("     Clinical Note: Patient is on Rivaroxaban, which can interfere with LA tests. However, the high and persistent values strongly suggest true LA positivity.")
    else:
        print("  -> LA laboratory criterion is NOT MET.")

    # 2b. Anticardiolipin (aCL) antibodies
    print("\n* Checking for Anticardiolipin (aCL) antibodies:")
    acl_igm_1_pos = labs1['aCL_IgM'] > criteria['norm_acl_ab2gp1']
    acl_igm_2_pos = labs2['aCL_IgM'] > criteria['norm_acl_ab2gp1']
    
    print(f"  - aCL IgM was {labs1['aCL_IgM']} UI/L 3 months ago and is now {labs2['aCL_IgM']} UI/L (Normal < {criteria['norm_acl_ab2gp1']} UI/L).")

    if acl_igm_1_pos and acl_igm_2_pos:
        print("  -> Finding: Persistent positive aCL IgM detected.")
        has_medium_titer = labs1['aCL_IgM'] > criteria['titer_acl_ab2gp1'] or labs2['aCL_IgM'] > criteria['titer_acl_ab2gp1']
        if has_medium_titer:
            print(f"     Titer is medium/high as one value ({labs2['aCL_IgM']}) is > {criteria['titer_acl_ab2gp1']} UI/L.")
            print("     -> aCL laboratory criterion is MET.")
            acl_met = True
        else:
            print("     -> aCL laboratory criterion is NOT MET as titer remains low.")
    else:
         print("     -> aCL laboratory criterion is NOT MET as it's not persistent.")
         
    # 2c. Anti-beta2-glycoprotein-I (aβ2GPI) antibodies
    print("\n* Checking for anti-ß2-glycoprotein-I (aß2GPI) antibodies:")
    ab2_igm_1_pos = labs1['ab2gp1_IgM'] > criteria['norm_acl_ab2gp1']
    ab2_igm_2_pos = labs2['ab2gp1_IgM'] > criteria['norm_acl_ab2gp1']

    print(f"  - aβ2GPI IgM was {labs1['ab2gp1_IgM']} UI/L 3 months ago and is now {labs2['ab2gp1_IgM']} UI/L (Normal < {criteria['norm_acl_ab2gp1']} UI/L).")

    if ab2_igm_1_pos and ab2_igm_2_pos:
        print("  -> Finding: Persistent positive aβ2GPI IgM detected.")
        has_high_titer = labs1['ab2gp1_IgM'] > criteria['titer_acl_ab2gp1'] or labs2['ab2gp1_IgM'] > criteria['titer_acl_ab2gp1']
        if has_high_titer:
             print(f"     Titer is >99th percentile as one value ({labs1['ab2gp1_IgM']}) is > {criteria['titer_acl_ab2gp1']} UI/L.")
             print("     -> aβ2GPI laboratory criterion is MET.")
             ab2gp1_met = True
        else:
             print("     -> aβ2GPI laboratory criterion is NOT MET as titer remains low.")
    else:
        print("     -> aβ2GPI laboratory criterion is NOT MET as it's not persistent.")
    
    lab_met = la_met or acl_met or ab2gp1_met
    print("-" * 30)

    # --- 3. Final Conclusion ---
    print("Step 3: Final Diagnosis")
    if clinical_met and lab_met:
        final_answer = "Yes"
        print("Conclusion: The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion.")
        print("This confirms the diagnosis of Antiphospholipid Syndrome (APS).")
        if la_met and acl_met and ab2gp1_met:
             print("Furthermore, the patient is 'triple positive' (positive for LA, aCL, and aß2GPI), which is associated with a very high risk for future thrombotic events.")
    else:
        final_answer = "No"
        print("Conclusion: The patient does not meet the full criteria for a definite diagnosis of Antiphospholipid Syndrome.")

    print(f"\nDoes this patient categorizes as having antiphospholipid syndrome?\n<<<{final_answer}>>>")

solve_aps_diagnosis()