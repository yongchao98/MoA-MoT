def check_antiphospholipid_syndrome():
    """
    Evaluates if a patient meets the criteria for Antiphospholipid Syndrome (APS)
    based on the revised Sapporo (Sydney) criteria.
    """
    # Patient Data
    patient_data = {
        "clinical": {
            "vte_events": 3,
            "pregnancy_morbidity": False # Not enough information for a positive criterion
        },
        "labs_3_months_ago": {
            "dRVVT_ratio": 1.44, "dRVVT_norm": 1.2,
            "aCL_IgM": 32, "aCL_norm": 20, "aCL_medium_titer": 40,
            "anti_beta2GPI_IgM": 41, "anti_beta2GPI_norm": 20
        },
        "labs_today": {
            "dRVVT_ratio": 1.51, "dRVVT_norm": 1.2,
            "aCL_IgM": 47, "aCL_norm": 20, "aCL_medium_titer": 40,
            "anti_beta2GPI_IgM": 29, "anti_beta2GPI_norm": 20
        },
        "lab_interval_weeks": 12
    }

    # --- Step 1: Check Clinical Criteria ---
    clinical_criterion_met = False
    print("Evaluating Clinical Criteria:")
    if patient_data["clinical"]["vte_events"] >= 1:
        clinical_criterion_met = True
        print(f"- Clinical criterion MET: Patient has {patient_data['clinical']['vte_events']} documented VTE events (>= 1 required).")
    else:
        print("- Clinical criterion NOT MET: No documented vascular thrombosis or relevant pregnancy morbidity.")

    # --- Step 2: Check Laboratory Criteria ---
    lab_criterion_met = False
    print("\nEvaluating Laboratory Criteria (must be positive on 2+ occasions >= 12 weeks apart):")
    
    labs1 = patient_data["labs_3_months_ago"]
    labs2 = patient_data["labs_today"]

    # Check for Lupus Anticoagulant (LA)
    la_positive_1 = labs1["dRVVT_ratio"] > labs1["dRVVT_norm"]
    la_positive_2 = labs2["dRVVT_ratio"] > labs2["dRVVT_norm"]
    persistent_la = la_positive_1 and la_positive_2
    if persistent_la:
        lab_criterion_met = True
        print(f"- Laboratory criterion MET (Lupus Anticoagulant): dRVVT ratio was {labs1['dRVVT_ratio']} and is now {labs2['dRVVT_ratio']} (Normal < {labs1['dRVVT_norm']}). This is persistently positive.")

    # Check for Anticardiolipin (aCL) antibodies
    acl_positive_1 = labs1["aCL_IgM"] > labs1["aCL_norm"]
    acl_positive_2 = labs2["aCL_IgM"] > labs2["aCL_norm"]
    # Check for medium/high titer
    medium_titer_present = labs1["aCL_IgM"] > labs1["aCL_medium_titer"] or labs2["aCL_IgM"] > labs2["aCL_medium_titer"]
    persistent_acl = acl_positive_1 and acl_positive_2 and medium_titer_present
    if persistent_acl:
        lab_criterion_met = True
        print(f"- Laboratory criterion MET (Anticardiolipin IgM): aCL IgM was {labs1['aCL_IgM']} and is now {labs2['aCL_IgM']} (Normal < {labs1['aCL_norm']}). This is persistently positive at a medium/high titer.")

    # Check for Anti-beta2-glycoprotein-I (anti-B2GPI) antibodies
    b2gpi_positive_1 = labs1["anti_beta2GPI_IgM"] > labs1["anti_beta2GPI_norm"]
    b2gpi_positive_2 = labs2["anti_beta2GPI_IgM"] > labs2["anti_beta2GPI_norm"]
    persistent_b2gpi = b2gpi_positive_1 and b2gpi_positive_2
    if persistent_b2gpi:
        lab_criterion_met = True
        print(f"- Laboratory criterion MET (Anti-B2GPI IgM): Anti-B2GPI IgM was {labs1['anti_beta2GPI_IgM']} and is now {labs2['anti_beta2GPI_IgM']} (Normal < {labs1['anti_beta2GPI_norm']}). This is persistently positive.")

    if not lab_criterion_met:
        print("- No persistent laboratory criterion was met.")

    # --- Step 3: Final Conclusion ---
    print("\n--- Conclusion ---")
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion AND at least one laboratory criterion.")
        print("Diagnosis of Antiphospholipid Syndrome is confirmed.")
    else:
        final_answer = "No"
        print("The patient does not meet the full criteria for Antiphospholipid Syndrome.")
    
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?\n<<<{final_answer}>>>")

if __name__ == "__main__":
    check_antiphospholipid_syndrome()