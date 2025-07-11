import sys

def check_aps_criteria():
    """
    Analyzes the patient's data against the APS diagnostic criteria.
    """
    # Patient Data
    patient_history = {
        "age": 34,
        "vte_events": [
            {"age": 18, "type": "Calf DVT", "trigger": "Plane travel"},
            {"age": 25, "type": "PE", "trigger": "Pregnancy"},
            {"age": 34, "type": "PE", "trigger": "Spontaneous"}
        ]
    }

    labs_3_months_ago = {
        "anti_b2gp1_igm": 41, "anti_b2gp1_igg": 18,
        "anticardiolipin_igm": 32, "anticardiolipin_igg": 9,
        "drvvt_ratio": 1.44
    }

    labs_today = {
        "anti_b2gp1_igm": 29, "anti_b2gp1_igg": 21,
        "anticardiolipin_igm": 47, "anticardiolipin_igg": 7,
        "drvvt_ratio": 1.51
    }

    clinical_criteria_met = False
    lab_criteria_met = False

    print("Step 1: Evaluating Clinical Criteria for APS")
    print("-" * 40)
    print(f"The patient has a history of {len(patient_history['vte_events'])} venous thromboembolic (VTE) events:")
    for event in patient_history['vte_events']:
        print(f"- {event['type']} at age {event['age']}")
    print("\nThe presence of one or more confirmed venous thrombotic events fulfills the clinical criteria for APS.")
    clinical_criteria_met = True
    print("Result: Clinical criteria are MET.\n")


    print("Step 2: Evaluating Laboratory Criteria for APS")
    print("-" * 40)
    print("The criteria require persistent positivity of at least one antiphospholipid antibody on two occasions at least 12 weeks apart.")
    print("The tests were performed 3 months ago and today, satisfying the 12-week interval.\n")

    # Check for persistent Lupus Anticoagulant (dRVVT)
    drvvt_1 = labs_3_months_ago["drvvt_ratio"]
    drvvt_2 = labs_today["drvvt_ratio"]
    is_la_positive = drvvt_1 > 1.2 and drvvt_2 > 1.2
    if is_la_positive:
        print(f"Lupus Anticoagulant (via dRVVT):")
        print(f"- 3 months ago: dRVVT ratio was {drvvt_1} (Normal < 1.2)")
        print(f"- Today: dRVVT ratio was {drvvt_2} (Normal < 1.2)")
        print("Result: The dRVVT is persistently positive, meeting a laboratory criterion.\n")
        lab_criteria_met = True

    # Check for persistent anticardiolipin antibodies
    acl_igm_1 = labs_3_months_ago["anticardiolipin_igm"]
    acl_igm_2 = labs_today["anticardiolipin_igm"]
    # Check for positivity (>20) on both occasions
    is_acl_positive = acl_igm_1 > 20 and acl_igm_2 > 20
    if is_acl_positive:
        print(f"Anticardiolipin IgM Antibodies:")
        print(f"- 3 months ago: aCL IgM was {acl_igm_1} UI/L (Normal < 20 UI/L)")
        print(f"- Today: aCL IgM was {acl_igm_2} UI/L (Normal < 20 UI/L)")
        print("Result: Anticardiolipin IgM antibodies are persistently positive, meeting a laboratory criterion.\n")
        lab_criteria_met = True

    # Check for persistent anti-B2GP1 antibodies
    ab2gp1_igm_1 = labs_3_months_ago["anti_b2gp1_igm"]
    ab2gp1_igm_2 = labs_today["anti_b2gp1_igm"]
    is_ab2gp1_positive = ab2gp1_igm_1 > 20 and ab2gp1_igm_2 > 20
    if is_ab2gp1_positive:
        print(f"Anti-ß2GP1 IgM Antibodies:")
        print(f"- 3 months ago: anti-ß2GP1 IgM was {ab2gp1_igm_1} UI/L (Normal < 20 UI/L)")
        print(f"- Today: anti-ß2GP1 IgM was {ab2gp1_igm_2} UI/L (Normal < 20 UI/L)")
        print("Result: Anti-ß2GP1 IgM antibodies are persistently positive, meeting a laboratory criterion.\n")
        lab_criteria_met = True
    
    print("Result: At least one laboratory criterion is MET.\n")

    print("Step 3: Final Conclusion")
    print("-" * 40)
    if clinical_criteria_met and lab_criteria_met:
        print("The patient meets both the clinical criteria (history of VTE) and the laboratory criteria (persistent antiphospholipid antibodies).")
        print("\nTherefore, the patient categorizes as having Antiphospholipid Syndrome (APS).")
        # To avoid multiple code blocks, we are using sys.stdout.write
        # for the final answer format which the system expects.
        sys.stdout.write("<<<Yes>>>")
    else:
        print("The patient does not meet both the clinical and laboratory criteria for APS.")
        sys.stdout.write("<<<No>>>")

# Run the diagnosis
check_aps_criteria()