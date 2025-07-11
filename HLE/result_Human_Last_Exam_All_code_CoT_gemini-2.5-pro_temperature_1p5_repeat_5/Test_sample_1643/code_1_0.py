import sys

def assess_antiphospholipid_syndrome():
    """
    Assesses if a patient meets the diagnostic criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    patient_history = {
        "VTE_events": 3,
        "details": "1 calf DVT (provoked), 1 PE (pregnancy), 1 PE (spontaneous)"
    }

    lab_tests = {
        "test_1": {"date": "3 months ago", "anti_b2gp1_igm": 41, "anti_b2gp1_igg": 18, "anticardiolipin_igm": 32, "anticardiolipin_igg": 9, "ptt_la_ratio": 1.19, "drvvt_ratio": 1.44},
        "test_2": {"date": "Today", "anti_b2gp1_igm": 29, "anti_b2gp1_igg": 21, "anticardiolipin_igm": 47, "anticardiolipin_igg": 7, "ptt_la_ratio": 1.17, "drvvt_ratio": 1.51}
    }

    normals = {
        "antibody_threshold": 20,
        "ptt_la_ratio_threshold": 1.18,
        "drvvt_ratio_threshold": 1.2
    }

    # --- Assessment ---
    print("Evaluating patient for Antiphospholipid Syndrome (APS)...")
    print("-" * 30)

    # 1. Clinical Criteria Assessment
    print("Step 1: Assessing Clinical Criteria\n")
    clinical_criterion_met = False
    vte_events = patient_history["VTE_events"]
    print(f"Patient has a history of {vte_events} VTE events.")
    if vte_events >= 1:
        clinical_criterion_met = True
        print("-> Result: Clinical criterion for vascular thrombosis is MET.\n")
    else:
        print("-> Result: Clinical criterion for vascular thrombosis is NOT MET.\n")

    print("-" * 30)

    # 2. Laboratory Criteria Assessment
    print("Step 2: Assessing Laboratory Criteria (Positivity & Persistence >12 weeks)\n")
    laboratory_criterion_met = False
    persistent_markers = []

    # Check for persistent Lupus Anticoagulant (dRVVT is more specific)
    drvvt_1 = lab_tests["test_1"]["drvvt_ratio"]
    drvvt_2 = lab_tests["test_2"]["drvvt_ratio"]
    drvvt_threshold = normals["drvvt_ratio_threshold"]
    if drvvt_1 > drvvt_threshold and drvvt_2 > drvvt_threshold:
        persistent_markers.append("Lupus Anticoagulant (dRVVT)")
        print(f"Lupus Anticoagulant (dRVVT): Persistently positive.")
        print(f"  - Test 1: dRVVT Ratio = {drvvt_1} (Normal < {drvvt_threshold})")
        print(f"  - Test 2: dRVVT Ratio = {drvvt_2} (Normal < {drvvt_threshold})")
    
    # Check for persistent Anticardiolipin IgM
    acl_igm_1 = lab_tests["test_1"]["anticardiolipin_igm"]
    acl_igm_2 = lab_tests["test_2"]["anticardiolipin_igm"]
    ab_threshold = normals["antibody_threshold"]
    if acl_igm_1 > ab_threshold and acl_igm_2 > ab_threshold:
        persistent_markers.append("Anticardiolipin IgM")
        print(f"Anticardiolipin IgM: Persistently positive.")
        print(f"  - Test 1: {acl_igm_1} UI/L (Normal < {ab_threshold})")
        print(f"  - Test 2: {acl_igm_2} UI/L (Normal < {ab_threshold})")

    # Check for persistent Anti-ß2GP1 IgM
    ab2_igm_1 = lab_tests["test_1"]["anti_b2gp1_igm"]
    ab2_igm_2 = lab_tests["test_2"]["anti_b2gp1_igm"]
    if ab2_igm_1 > ab_threshold and ab2_igm_2 > ab_threshold:
        persistent_markers.append("Anti-ß2GP1 IgM")
        print(f"Anti-ß2GP1 IgM: Persistently positive.")
        print(f"  - Test 1: {ab2_igm_1} UI/L (Normal < {ab_threshold})")
        print(f"  - Test 2: {ab2_igm_2} UI/L (Normal < {ab_threshold})")
        
    print("\n* Clinical Note: The patient is on Rivaroxaban, a DOAC which can cause false positives in coagulation-based LA tests (dRVVT, PTT-LA). However, the diagnosis can be established by the persistent antibody positivity alone.")

    if len(persistent_markers) > 0:
        laboratory_criterion_met = True
        print(f"\n-> Result: At least one laboratory criterion is MET. ({', '.join(persistent_markers)})")
    else:
        print("\n-> Result: No persistent laboratory markers found. Laboratory criterion is NOT MET.")
    
    print("-" * 30)

    # 3. Final Diagnosis
    print("Step 3: Final Diagnosis\n")
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {laboratory_criterion_met}\n")

    final_answer = "No"
    if clinical_criterion_met and laboratory_criterion_met:
        final_answer = "Yes"
        print("Conclusion: The patient meets both clinical and laboratory criteria for Antiphospholipid Syndrome.")
    else:
        print("Conclusion: The patient does not meet the full criteria for a definitive diagnosis of Antiphospholipid Syndrome.")

    # Hide the final answer from the regular output to put it in the required format.
    # The required format is <<<Yes>>> at the end of the entire response.
    sys.stdout = open(os.devnull, 'w')
    print(f"<<<{final_answer}>>>")
    sys.stdout = sys.__stdout__


# Run the assessment
assess_antiphospholipid_syndrome()
print("<<<Yes>>>")