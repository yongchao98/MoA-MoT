def check_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) and prints the reasoning.
    """

    # Patient Data
    patient = {
        "clinical": {
            "vte_events": 3,
            "pregnancy_morbidity": False # Not enough information for a positive criterion
        },
        "labs_3_months_ago": {
            "aCL_IgM": 32, "aCL_IgG": 9,
            "ab2gp1_IgM": 41, "ab2gp1_IgG": 18,
            "ptt_la_ratio": 1.19,
            "drvvT_ratio": 1.44
        },
        "labs_today": {
            "aCL_IgM": 47, "aCL_IgG": 7,
            "ab2gp1_IgM": 29, "ab2gp1_IgG": 21,
            "ptt_la_ratio": 1.17,
            "drvvT_ratio": 1.51
        },
        "normals": {
            "aCL_IgM": 20, "aCL_IgG": 20,
            "ab2gp1_IgM": 20, "ab2gp1_IgG": 20,
            "ptt_la_ratio": 1.18,
            "drvvT_ratio": 1.2
        }
    }

    print("Evaluating patient for Antiphospholipid Syndrome (APS)...")
    print("-" * 50)

    # Step 1: Check Clinical Criteria
    print("Step 1: Checking Clinical Criteria (requires >= 1 vascular thrombosis)")
    clinical_criterion_met = False
    if patient["clinical"]["vte_events"] >= 1:
        clinical_criterion_met = True
        print(f"Result: Clinical criterion MET. Patient has a history of {patient['clinical']['vte_events']} VTE events.")
    else:
        print("Result: Clinical criterion NOT MET.")

    print("-" * 50)

    # Step 2: Check Laboratory Criteria
    print("Step 2: Checking Laboratory Criteria (requires >= 1 persistent positive antibody test >12 weeks apart)")
    lab_criterion_met = False

    # Check for persistent Lupus Anticoagulant (LA)
    # Note: DOACs like Rivaroxaban can cause false positives, but dRVVT is strongly positive.
    la_test1_positive = patient["labs_3_months_ago"]["drvvT_ratio"] > patient["normals"]["drvvT_ratio"]
    la_test2_positive = patient["labs_today"]["drvvT_ratio"] > patient["normals"]["drvvT_ratio"]
    if la_test1_positive and la_test2_positive:
        lab_criterion_met = True
        print(f"-> Lupus Anticoagulant: POSITIVE")
        print(f"   - 3 months ago: dRVVT ratio was {patient['labs_3_months_ago']['drvvT_ratio']} (Normal < {patient['normals']['drvvT_ratio']})")
        print(f"   - Today: dRVVT ratio is {patient['labs_today']['drvvT_ratio']} (Normal < {patient['normals']['drvvT_ratio']})")
        print("   This is persistently positive.")
    else:
        print("-> Lupus Anticoagulant: Negative or not persistent.")

    # Check for persistent Anticardiolipin (aCL) antibodies
    acl_test1_positive = patient["labs_3_months_ago"]["aCL_IgM"] > patient["normals"]["aCL_IgM"] or \
                         patient["labs_3_months_ago"]["aCL_IgG"] > patient["normals"]["aCL_IgG"]
    acl_test2_positive = patient["labs_today"]["aCL_IgM"] > patient["normals"]["aCL_IgM"] or \
                         patient["labs_today"]["aCL_IgG"] > patient["normals"]["aCL_IgG"]
    if acl_test1_positive and acl_test2_positive:
        lab_criterion_met = True
        print(f"-> Anticardiolipin Antibodies: POSITIVE")
        print(f"   - 3 months ago: aCL IgM was {patient['labs_3_months_ago']['aCL_IgM']} (Normal < {patient['normals']['aCL_IgM']})")
        print(f"   - Today: aCL IgM is {patient['labs_today']['aCL_IgM']} (Normal < {patient['normals']['aCL_IgM']})")
        print("   This is persistently positive.")
    else:
        print("-> Anticardiolipin Antibodies: Negative or not persistent.")

    # Check for persistent Anti-B2GP1 antibodies
    ab2gp1_test1_positive = patient["labs_3_months_ago"]["ab2gp1_IgM"] > patient["normals"]["ab2gp1_IgM"] or \
                            patient["labs_3_months_ago"]["ab2gp1_IgG"] > patient["normals"]["ab2gp1_IgG"]
    ab2gp1_test2_positive = patient["labs_today"]["ab2gp1_IgM"] > patient["normals"]["ab2gp1_IgM"] or \
                            patient["labs_today"]["ab2gp1_IgG"] > patient["normals"]["ab2gp1_IgG"]
    if ab2gp1_test1_positive and ab2gp1_test2_positive:
        lab_criterion_met = True
        print(f"-> Anti-B2GP1 Antibodies: POSITIVE")
        print(f"   - 3 months ago: anti-B2GP1 IgM was {patient['labs_3_months_ago']['ab2gp1_IgM']} (Normal < {patient['normals']['ab2gp1_IgM']})")
        print(f"   - Today: anti-B2GP1 IgM/IgG are {patient['labs_today']['ab2gp1_IgM']}/{patient['labs_today']['ab2gp1_IgG']} respectively (Normal < {patient['normals']['ab2gp1_IgM']})")
        print("   This is persistently positive.")
    else:
        print("-> Anti-B2GP1 Antibodies: Negative or not persistent.")

    if lab_criterion_met:
        print("\nResult: Laboratory criterion MET. At least one antibody is persistently positive.")
    else:
        print("\nResult: Laboratory criterion NOT MET.")

    print("-" * 50)

    # Step 3: Final Diagnosis
    print("Step 3: Final Conclusion")
    final_answer = "No"
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets both clinical and laboratory criteria for APS.")
    else:
        print("The patient does not meet both clinical and laboratory criteria for APS.")

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome? \nAnswer: {final_answer}")
    return final_answer

# Run the analysis and print the final answer in the required format
final_diagnosis = check_aps_diagnosis()
print(f"<<<{final_diagnosis}>>>")