def check_antiphospholipid_syndrome():
    """
    Analyzes patient data against the 2006 Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    patient_data = {
        "history": [
            "Calf DVT at age 18",
            "Intermediate-risk PE at age 25",
            "PE 4 months ago, spontaneous"
        ],
        "labs_3_months_ago": {
            "anti_b2GP1_IgM": 41,  # N < 20 UI/L
            "anti_b2GP1_IgG": 18,  # N < 20 UI/L
            "anticardiolipin_IgM": 32, # N < 20 UI/L
            "anticardiolipin_IgG": 9,  # N < 20 UI/L
            "dRVVT_ratio": 1.44     # N < 1.2
        },
        "labs_today": {
            "anti_b2GP1_IgM": 29,  # N < 20 UI/L
            "anti_b2GP1_IgG": 21,  # N < 20 UI/L
            "anticardiolipin_IgM": 47, # N < 20 UI/L
            "anticardiolipin_IgG": 7,  # N < 20 UI/L
            "dRVVT_ratio": 1.51     # N < 1.2
        }
    }

    # --- Step 1: Evaluate Clinical Criteria ---
    print("--- Evaluating Clinical Criteria for APS ---")
    thrombotic_events = len(patient_data["history"])
    clinical_criterion_met = thrombotic_events >= 1

    print(f"The patient has a history of {thrombotic_events} venous thromboembolic (VTE) events:")
    for event in patient_data["history"]:
        print(f"- {event}")
    
    if clinical_criterion_met:
        print("\nConclusion: The patient meets the clinical criterion for APS (Vascular Thrombosis).\n")
    else:
        print("\nConclusion: The patient does not meet the clinical criterion for APS.\n")

    # --- Step 2: Evaluate Laboratory Criteria ---
    print("--- Evaluating Laboratory Criteria for APS ---")
    print("Criteria requires at least one positive lab test, confirmed on two occasions at least 12 weeks apart.")
    print("The tests were performed 3 months ago and today, which satisfies the 12-week interval.\n")
    
    labs1 = patient_data["labs_3_months_ago"]
    labs2 = patient_data["labs_today"]
    lab_criteria_met = []

    # Check for Lupus Anticoagulant (LA) via dRVVT
    la_positive_1 = labs1["dRVVT_ratio"] > 1.2
    la_positive_2 = labs2["dRVVT_ratio"] > 1.2
    if la_positive_1 and la_positive_2:
        la_criterion_met = True
        lab_criteria_met.append("Lupus Anticoagulant")
        print(f"1. Lupus Anticoagulant (LA): POSITIVE")
        print(f"   - 3 months ago: dRVVT ratio was {labs1['dRVVT_ratio']} (Positive)")
        print(f"   - Today: dRVVT ratio is {labs2['dRVVT_ratio']} (Positive)")
        print("   - Result: Criterion MET.\n")
    else:
        la_criterion_met = False
        print("1. Lupus Anticoagulant (LA): NEGATIVE\n")


    # Check for Anticardiolipin (aCL) antibodies (IgG or IgM)
    # Criterion requires persistence and medium-to-high titer (>40 MPL/GPL or >99th percentile)
    acl_igm_persistent = labs1["anticardiolipin_IgM"] > 20 and labs2["anticardiolipin_IgM"] > 20
    acl_igg_persistent = labs1["anticardiolipin_IgG"] > 20 and labs2["anticardiolipin_IgG"] > 20
    medium_high_titer = labs1["anticardiolipin_IgM"] > 40 or labs2["anticardiolipin_IgM"] > 40 or \
                        labs1["anticardiolipin_IgG"] > 40 or labs2["anticardiolipin_IgG"] > 40

    if (acl_igm_persistent or acl_igg_persistent) and medium_high_titer:
        acl_criterion_met = True
        lab_criteria_met.append("Anticardiolipin Antibodies")
        print(f"2. Anticardiolipin Antibodies (aCL): POSITIVE")
        print(f"   - Anticardiolipin IgM was persistently positive ({labs1['anticardiolipin_IgM']} and {labs2['anticardiolipin_IgM']} U/L).")
        print(f"   - The titer today ({labs2['anticardiolipin_IgM']} U/L) is medium/high (>40).")
        print("   - Result: Criterion MET.\n")
    else:
        acl_criterion_met = False
        print("2. Anticardiolipin Antibodies (aCL): NEGATIVE\n")

    # Check for Anti-β2 Glycoprotein I (aβ2GPI) antibodies (IgG or IgM)
    ab2gpi_igm_persistent = labs1["anti_b2GP1_IgM"] > 20 and labs2["anti_b2GP1_IgM"] > 20
    ab2gpi_igg_persistent = labs1["anti_b2GP1_IgG"] > 20 and labs2["anti_b2GP1_IgG"] > 20

    if ab2gpi_igm_persistent or ab2gpi_igg_persistent:
        ab2gpi_criterion_met = True
        lab_criteria_met.append("Anti-β2 Glycoprotein I Antibodies")
        print(f"3. Anti-β2 Glycoprotein I Antibodies (aβ2GPI): POSITIVE")
        print(f"   - Anti-β2GPI IgM was persistently positive ({labs1['anti_b2GP1_IgM']} and {labs2['anti_b2GP1_IgM']} U/L).")
        print(f"   - Anti-β2GPI IgG became positive today ({labs2['anti_b2GP1_IgG']} U/L).")
        print("   - Result: Criterion MET.\n")
    else:
        ab2gpi_criterion_met = False
        print("3. Anti-β2 Glycoprotein I Antibodies (aβ2GPI): NEGATIVE\n")

    any_lab_criterion_met = la_criterion_met or acl_criterion_met or ab2gpi_criterion_met
    if any_lab_criterion_met:
        print(f"Conclusion: The patient meets {len(lab_criteria_met)} laboratory criteria for APS ({', '.join(lab_criteria_met)}).\n")
    else:
        print("Conclusion: The patient does not meet any laboratory criteria for APS.\n")

    # --- Step 3: Final Diagnosis ---
    print("--- Final Diagnosis ---")
    if clinical_criterion_met and any_lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion (in this case, all three).")
        print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    else:
        final_answer = "No"
        print("The patient does not meet both the minimum clinical and laboratory criteria for APS.")
        print("\nDoes this patient categorize as having antiphospholipid syndrome?")

    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    check_antiphospholipid_syndrome()