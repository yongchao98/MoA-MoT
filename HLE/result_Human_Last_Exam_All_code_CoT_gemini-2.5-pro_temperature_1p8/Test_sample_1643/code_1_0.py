def check_aps_diagnosis():
    """
    Evaluates if a patient meets the criteria for Antiphospholipid Syndrome (APS)
    based on their clinical history and laboratory results.
    """

    # --- Patient Data ---
    clinical_history = {
        "vte_events": 3
    }

    lab_3_months_ago = {
        "aCL_IgM": 32, "aCL_IgG": 9,
        "aB2GPI_IgM": 41, "aB2GPI_IgG": 18,
        "dRVVT_ratio": 1.44
    }

    lab_today = {
        "aCL_IgM": 47, "aCL_IgG": 7,
        "aB2GPI_IgM": 29, "aB2GPI_IgG": 21,
        "dRVVT_ratio": 1.51
    }

    # Normal thresholds
    # Note: For aCL and aB2GPI, the criteria mention "medium or high titer".
    # Given the lab normal of <20, any value >20 is considered positive for this evaluation.
    # For dRVVT, N < 1.2
    lab_normals = {
        "antibody_threshold": 20,
        "dRVVT_threshold": 1.2
    }

    print("Step 1: Evaluating Clinical Criteria for APS...")
    # --- Check Clinical Criteria ---
    clinical_criterion_met = False
    if clinical_history["vte_events"] >= 1:
        clinical_criterion_met = True
        print(f"Outcome: MET. The patient has had {clinical_history['vte_events']} documented VTE event(s). This satisfies the clinical criterion of vascular thrombosis.\n")
    else:
        print("Outcome: NOT MET. The patient has no documented VTE events.\n")

    print("Step 2: Evaluating Laboratory Criteria for APS (Persistence over 12 weeks)...")
    # --- Check Laboratory Criteria ---
    la_positive = False
    acl_positive = False
    ab2gpi_positive = False

    # Check for Lupus Anticoagulant (LA) persistence using dRVVT
    print("Checking Lupus Anticoagulant (dRVVT)...")
    val1 = lab_3_months_ago['dRVVT_ratio']
    val2 = lab_today['dRVVT_ratio']
    threshold = lab_normals['dRVVT_threshold']
    if val1 > threshold and val2 > threshold:
        la_positive = True
        print(f"Result: POSITIVE. dRVVT ratio was persistently elevated ({val1} and {val2}, both > {threshold}).")
    else:
        print(f"Result: NEGATIVE. dRVVT ratio was not persistently elevated ({val1} and {val2} vs. threshold {threshold}).")

    # Check for anticardiolipin (aCL) persistence
    print("\nChecking anticardiolipin antibodies (aCL)...")
    acl_igm_persistent = lab_3_months_ago['aCL_IgM'] > lab_normals['antibody_threshold'] and lab_today['aCL_IgM'] > lab_normals['antibody_threshold']
    acl_igg_persistent = lab_3_months_ago['aCL_IgG'] > lab_normals['antibody_threshold'] and lab_today['aCL_IgG'] > lab_normals['antibody_threshold']
    
    if acl_igm_persistent or acl_igg_persistent:
        acl_positive = True
        if acl_igm_persistent:
            print(f"Result: POSITIVE. aCL IgM was persistently positive ({lab_3_months_ago['aCL_IgM']} and {lab_today['aCL_IgM']}, both > {lab_normals['antibody_threshold']}).")
        if acl_igg_persistent:
             print(f"Result: POSITIVE. aCL IgG was persistently positive ({lab_3_months_ago['aCL_IgG']} and {lab_today['aCL_IgG']}, both > {lab_normals['antibody_threshold']}).")
    else:
        print("Result: NEGATIVE. No persistent aCL antibodies detected.")

    # Check for anti-β2-glycoprotein I (aβ2GPI) persistence
    print("\nChecking anti-β2-glycoprotein I antibodies (aβ2GPI)...")
    ab2gpi_igm_persistent = lab_3_months_ago['aB2GPI_IgM'] > lab_normals['antibody_threshold'] and lab_today['aB2GPI_IgM'] > lab_normals['antibody_threshold']
    ab2gpi_igg_persistent = lab_3_months_ago['aB2GPI_IgG'] > lab_normals['antibody_threshold'] and lab_today['aB2GPI_IgG'] > lab_normals['antibody_threshold']

    if ab2gpi_igm_persistent or ab2gpi_igg_persistent:
        ab2gpi_positive = True
        if ab2gpi_igm_persistent:
            print(f"Result: POSITIVE. aβ2GPI IgM was persistently positive ({lab_3_months_ago['aB2GPI_IgM']} and {lab_today['aB2GPI_IgM']}, both > {lab_normals['antibody_threshold']}).")
        if ab2gpi_igg_persistent:
             print(f"Result: POSITIVE. aβ2GPI IgG was persistently positive ({lab_3_months_ago['aB2GPI_IgG']} and {lab_today['aB2GPI_IgG']}, both > {lab_normals['antibody_threshold']}).")
    else:
        print("Result: NEGATIVE. No persistent aβ2GPI antibodies detected.")

    laboratory_criterion_met = la_positive or acl_positive or ab2gpi_positive
    
    if laboratory_criterion_met:
        print("\nOutcome: MET. At least one laboratory criterion was met.\n")
    else:
        print("\nOutcome: NOT MET. No laboratory criteria were met.\n")

    print("Step 3: Final Conclusion...")
    if clinical_criterion_met and laboratory_criterion_met:
        final_answer = "Yes"
        print("Diagnosis: The patient meets both clinical (vascular thrombosis) and laboratory criteria.")
    else:
        final_answer = "No"
        print("Diagnosis: The patient does not meet the full criteria for APS.")
    
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

check_aps_diagnosis()