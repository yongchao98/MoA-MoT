def analyze_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the 2006 Sydney consensus criteria.
    """

    # --- Patient Data ---
    clinical_events = [
        {"type": "Calf DVT", "age": 18},
        {"type": "PE", "age": 25},
        {"type": "PE", "age": 34}
    ]

    lab_3_months_ago = {
        "aCL_IgM": 32, "aCL_IgG": 9,
        "antiB2GPI_IgM": 41, "antiB2GPI_IgG": 18,
        "dRVVT_ratio": 1.44
    }

    lab_today = {
        "aCL_IgM": 47, "aCL_IgG": 7,
        "antiB2GPI_IgM": 29, "antiB2GPI_IgG": 21,
        "dRVVT_ratio": 1.51
    }

    print("Evaluating patient for Antiphospholipid Syndrome (APS)...")
    print("Diagnosis requires at least ONE clinical criterion AND ONE laboratory criterion.\n")

    # --- Step 1: Evaluate Clinical Criteria ---
    print("--- Step 1: Clinical Criteria Evaluation ---")
    print("Checking for a history of one or more vascular thrombotic events.")
    
    thrombotic_event_count = len(clinical_events)
    clinical_criterion_met = thrombotic_event_count >= 1
    
    if clinical_criterion_met:
        print(f"Result: MET. Patient has a history of {thrombotic_event_count} VTE events.")
    else:
        print("Result: NOT MET.")
    
    # --- Step 2: Evaluate Laboratory Criteria ---
    print("\n--- Step 2: Laboratory Criteria Evaluation ---")
    print("Checking for persistent antiphospholipid antibodies (tests were 3 months apart, satisfying >=12 week interval).\n")
    
    lab_criteria_met = []

    # Check for persistent Lupus Anticoagulant (LA) via dRVVT
    dRVVT_1, dRVVT_2 = lab_3_months_ago["dRVVT_ratio"], lab_today["dRVVT_ratio"]
    if dRVVT_1 > 1.2 and dRVVT_2 > 1.2:
        lab_criteria_met.append(True)
        print(f"Lupus Anticoagulant (dRVVT): PERSISTENTLY POSITIVE")
        print(f" -> dRVVT ratio was {dRVVT_1} and is now {dRVVT_2} (Normal < 1.2).")
    else:
        lab_criteria_met.append(False)
        print("Lupus Anticoagulant (dRVVT): Not persistently positive.")

    # Check for persistent anticardiolipin (aCL) antibodies (IgM)
    aCL_IgM_1, aCL_IgM_2 = lab_3_months_ago["aCL_IgM"], lab_today["aCL_IgM"]
    # Check for persistent positive (>20) AND at least one result is medium/high titer (>40)
    if aCL_IgM_1 > 20 and aCL_IgM_2 > 20 and (aCL_IgM_1 > 40 or aCL_IgM_2 > 40):
        lab_criteria_met.append(True)
        print(f"Anticardiolipin IgM (aCL IgM): PERSISTENTLY POSITIVE (medium titer)")
        print(f" -> Value was {aCL_IgM_1} and is now {aCL_IgM_2} UI/L (Normal < 20).")
    else:
        lab_criteria_met.append(False)
        print(f"Anticardiolipin IgM (aCL IgM): Not meeting criteria for persistent medium/high titer.")

    # Check for persistent anti-B2GP1 antibodies (IgM)
    b2gp1_IgM_1, b2gp1_IgM_2 = lab_3_months_ago["antiB2GPI_IgM"], lab_today["antiB2GPI_IgM"]
    if b2gp1_IgM_1 > 20 and b2gp1_IgM_2 > 20:
        lab_criteria_met.append(True)
        print(f"Anti-B2GP1 IgM: PERSISTENTLY POSITIVE")
        print(f" -> Value was {b2gp1_IgM_1} and is now {b2gp1_IgM_2} UI/L (Normal < 20).")
    else:
        lab_criteria_met.append(False)
        print("Anti-B2GP1 IgM: Not persistently positive.")

    overall_lab_criterion_met = any(lab_criteria_met)
    print(f"\nResult: {'MET' if overall_lab_criterion_met else 'NOT MET'}. At least one laboratory criterion was met.")

    # --- Step 3: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    if clinical_criterion_met and overall_lab_criterion_met:
        print("The patient meets both clinical and laboratory criteria for a definitive diagnosis of Antiphospholipid Syndrome.")
        final_answer = "Yes"
    else:
        print("The patient does not meet the full criteria for a definitive diagnosis of Antiphospholipid Syndrome.")
        final_answer = "No"

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?\n<<<" + final_answer + ">>>")

if __name__ == "__main__":
    analyze_aps_diagnosis()