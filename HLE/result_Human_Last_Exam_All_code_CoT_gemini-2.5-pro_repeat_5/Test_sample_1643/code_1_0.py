def diagnose_aps():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the Sydney criteria.
    """

    # Patient Data
    patient = {
        "age": 34,
        "vte_events": [
            {"type": "Calf DVT", "age": 18, "trigger": "Plane travel"},
            {"type": "PE", "age": 25, "trigger": "Pregnancy"},
            {"type": "PE", "trigger": "Spontaneous"}
        ],
        "lab_tests": {
            "3_months_ago": {
                "aCL_IgM": 41, "aCL_IgG": 9,
                "anti_b2GP1_IgM": 32, "anti_b2GP1_IgG": 18,
                "dRVVT_ratio": 1.44
            },
            "today": {
                "aCL_IgM": 47, "aCL_IgG": 7,
                "anti_b2GP1_IgM": 29, "anti_b2GP1_IgG": 21,
                "dRVVT_ratio": 1.51
            }
        },
        "lab_normals": {
            "aCL_IgM": 20, "aCL_IgG": 20,
            "anti_b2GP1_IgM": 20, "anti_b2GP1_IgG": 20,
            "dRVVT_ratio": 1.2
        }
    }

    print("Step 1: Evaluating Clinical Criteria for APS")
    print("---------------------------------------------")

    # Check for Vascular Thrombosis
    num_vte_events = len(patient["vte_events"])
    clinical_criteria_met = num_vte_events >= 1

    print(f"Patient has experienced {num_vte_events} venous thromboembolic (VTE) events.")
    if clinical_criteria_met:
        print("Finding: The patient meets the clinical criterion of having one or more episodes of venous thrombosis.")
    else:
        print("Finding: The patient does not meet the clinical criterion for vascular thrombosis.")
    print("\n")


    print("Step 2: Evaluating Laboratory Criteria for APS")
    print("--------------------------------------------")
    print("Requirement: Persistence of at least one marker on two occasions >= 12 weeks apart.")
    print("The tests were performed 3 months ago and today, satisfying the time interval.\n")

    lab_criteria_met = False
    
    # Check for Lupus Anticoagulant (dRVVT)
    drvvt_1 = patient["lab_tests"]["3_months_ago"]["dRVVT_ratio"]
    drvvt_2 = patient["lab_tests"]["today"]["dRVVT_ratio"]
    drvvt_norm = patient["lab_normals"]["dRVVT_ratio"]
    la_positive = drvvt_1 > drvvt_norm and drvvt_2 > drvvt_norm
    if la_positive:
        print(f"1. Lupus Anticoagulant (LA):")
        print(f"   - 3 months ago: dRVVT ratio was {drvvt_1} (Normal < {drvvt_norm}) -> Positive")
        print(f"   - Today: dRVVT ratio is {drvvt_2} (Normal < {drvvt_norm}) -> Positive")
        print(f"   Finding: LA is persistently positive.")
        lab_criteria_met = True
    print("")

    # Check for Anticardiolipin Antibodies (aCL)
    # We check for medium/high titer, often defined as >40 or >99th percentile.
    # Here, we check for any persistent elevation above normal.
    acl_igm_1 = patient["lab_tests"]["3_months_ago"]["aCL_IgM"]
    acl_igm_2 = patient["lab_tests"]["today"]["aCL_IgM"]
    acl_norm = patient["lab_normals"]["aCL_IgM"]
    acl_igm_positive = acl_igm_1 > acl_norm and acl_igm_2 > acl_norm
    if acl_igm_positive:
        print(f"2. Anticardiolipin (aCL) IgM Antibodies:")
        print(f"   - 3 months ago: aCL IgM was {acl_igm_1} UI/L (Normal < {acl_norm}) -> Positive")
        print(f"   - Today: aCL IgM is {acl_igm_2} UI/L (Normal < {acl_norm}) -> Positive")
        print(f"   Finding: aCL IgM antibodies are persistently positive.")
        lab_criteria_met = True
    print("")
    
    # Check for Anti-B2GP1 Antibodies
    b2gp1_igm_1 = patient["lab_tests"]["3_months_ago"]["anti_b2GP1_IgM"]
    b2gp1_igm_2 = patient["lab_tests"]["today"]["anti_b2GP1_IgM"]
    b2gp1_norm = patient["lab_normals"]["anti_b2GP1_IgM"]
    b2gp1_igm_positive = b2gp1_igm_1 > b2gp1_norm and b2gp1_igm_2 > b2gp1_norm
    if b2gp1_igm_positive:
        print(f"3. Anti-B2-Glycoprotein-I (anti-B2GPI) IgM Antibodies:")
        print(f"   - 3 months ago: anti-B2GPI IgM was {b2gp1_igm_1} UI/L (Normal < {b2gp1_norm}) -> Positive")
        print(f"   - Today: anti-B2GPI IgM is {b2gp1_igm_2} UI/L (Normal < {b2gp1_norm}) -> Positive")
        print(f"   Finding: Anti-B2GPI IgM antibodies are persistently positive.")
        lab_criteria_met = True
    print("\n")


    print("Step 3: Final Conclusion")
    print("-------------------------")
    final_diagnosis = "No"
    if clinical_criteria_met and lab_criteria_met:
        final_diagnosis = "Yes"
        print("The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion (persistently positive antiphospholipid antibodies).")
        print("Therefore, the patient is categorized as having Antiphospholipid Syndrome.")
    elif clinical_criteria_met:
        print("The patient meets the clinical criteria, but the laboratory criteria for persistence are not met.")
    elif lab_criteria_met:
        print("The patient meets the laboratory criteria, but does not have a documented clinical event qualifying for APS.")
    else:
        print("The patient does not meet the clinical or laboratory criteria for APS.")

    return final_diagnosis

if __name__ == '__main__':
    result = diagnose_aps()
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?\n<<<{result}>>>")