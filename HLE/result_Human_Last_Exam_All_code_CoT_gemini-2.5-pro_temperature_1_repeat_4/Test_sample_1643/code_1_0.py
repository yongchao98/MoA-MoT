def diagnose_antiphospholipid_syndrome():
    """
    Analyzes patient data against the 2006 revised Sapporo criteria for APS.
    """

    # --- Patient Data ---
    patient_info = {
        "vte_events": 3,
    }

    lab_test_1 = {
        "time": "3 months ago",
        "anti_b2gp1_igm": 41, "anti_b2gp1_igm_norm": 20,
        "anticardiolipin_igm": 32, "anticardiolipin_igm_norm": 20,
        "drvvT_ratio": 1.44, "drvvT_ratio_norm": 1.2,
    }

    lab_test_2 = {
        "time": "Today",
        "anti_b2gp1_igm": 29, "anti_b2gp1_igm_norm": 20,
        "anticardiolipin_igm": 47, "anticardiolipin_igm_norm": 20,
        "drvvT_ratio": 1.51, "drvvT_ratio_norm": 1.2,
    }

    print("Evaluating patient for Antiphospholipid Syndrome (APS) based on the 2006 revised Sapporo criteria.")
    print("-" * 60)

    # --- Step 1: Assess Clinical Criteria ---
    print("Step 1: Assessing Clinical Criteria")
    clinical_criteria_met = patient_info["vte_events"] >= 1
    if clinical_criteria_met:
        print(f"-> Clinical Criterion MET: Patient has a history of {patient_info['vte_events']} vascular thrombosis events.")
    else:
        print("-> Clinical Criterion NOT MET.")
    print("-" * 60)


    # --- Step 2: Assess Laboratory Criteria (Persistence over >=12 weeks) ---
    print("Step 2: Assessing Laboratory Criteria")
    print("The tests are 3 months apart, satisfying the >=12 week interval requirement.\n")

    # Check for Lupus Anticoagulant (LA) via dRVVT
    la_t1_positive = lab_test_1["drvvT_ratio"] > lab_test_1["drvvT_ratio_norm"]
    la_t2_positive = lab_test_2["drvvT_ratio"] > lab_test_2["drvvT_ratio_norm"]
    persistent_la = la_t1_positive and la_t2_positive
    print("Assessing Lupus Anticoagulant (LA) via dRVVT ratio:")
    print(f"  - Test 1 (3 months ago): dRVVT = {lab_test_1['drvvT_ratio']} (Normal < {lab_test_1['drvvT_ratio_norm']}) -> {'Positive' if la_t1_positive else 'Negative'}")
    print(f"  - Test 2 (Today): dRVVT = {lab_test_2['drvvT_ratio']} (Normal < {lab_test_2['drvvT_ratio_norm']}) -> {'Positive' if la_t2_positive else 'Negative'}")
    if persistent_la:
        print("  -> Result: LA is persistently positive. (NOTE: This result is highly likely to be a false positive due to Rivaroxaban interference).")
    else:
        print("  -> Result: LA is not persistently positive.")
    print("")

    # Check for anticardiolipin (aCL) antibodies
    acl_t1_positive = lab_test_1["anticardiolipin_igm"] > lab_test_1["anticardiolipin_igm_norm"]
    acl_t2_positive = lab_test_2["anticardiolipin_igm"] > lab_test_2["anticardiolipin_igm_norm"]
    persistent_acl = acl_t1_positive and acl_t2_positive
    print("Assessing anticardiolipin (aCL) IgM antibodies:")
    print(f"  - Test 1 (3 months ago): aCL IgM = {lab_test_1['anticardiolipin_igm']} UI/L (Normal < {lab_test_1['anticardiolipin_igm_norm']}) -> {'Positive' if acl_t1_positive else 'Negative'}")
    print(f"  - Test 2 (Today): aCL IgM = {lab_test_2['anticardiolipin_igm']} UI/L (Normal < {lab_test_2['anticardiolipin_igm_norm']}) -> {'Positive' if acl_t2_positive else 'Negative'}")
    if persistent_acl:
        print("  -> Result: aCL IgM is persistently positive.")
    else:
        print("  -> Result: aCL IgM is not persistently positive.")
    print("")

    # Check for anti-ß2-glycoprotein-I (anti-ß2GPI) antibodies
    ab2_t1_positive = lab_test_1["anti_b2gp1_igm"] > lab_test_1["anti_b2gp1_igm_norm"]
    ab2_t2_positive = lab_test_2["anti_b2gp1_igm"] > lab_test_2["anti_b2gp1_igm_norm"]
    persistent_ab2 = ab2_t1_positive and ab2_t2_positive
    print("Assessing anti-ß2-glycoprotein-I (anti-ß2GPI) IgM antibodies:")
    print(f"  - Test 1 (3 months ago): anti-ß2GPI IgM = {lab_test_1['anti_b2gp1_igm']} UI/L (Normal < {lab_test_1['anti_b2gp1_igm_norm']}) -> {'Positive' if ab2_t1_positive else 'Negative'}")
    print(f"  - Test 2 (Today): anti-ß2GPI IgM = {lab_test_2['anti_b2gp1_igm']} UI/L (Normal < {lab_test_2['anti_b2gp1_igm_norm']}) -> {'Positive' if ab2_t2_positive else 'Negative'}")
    if persistent_ab2:
        print("  -> Result: anti-ß2GPI IgM is persistently positive.")
    else:
        print("  -> Result: anti-ß2GPI IgM is not persistently positive.")
    print("-" * 60)

    # --- Step 3: Final Conclusion ---
    laboratory_criteria_met = persistent_la or persistent_acl or persistent_ab2
    print("Step 3: Final Conclusion")
    if clinical_criteria_met and laboratory_criteria_met:
        print("-> The patient meets at least one clinical criterion (Vascular Thrombosis) and at least one laboratory criterion (persistently positive aCL and anti-ß2GPI antibodies).")
        final_answer = "Yes"
    else:
        print("-> The patient does not meet both the clinical and laboratory criteria.")
        final_answer = "No"

    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    diagnose_antiphospholipid_syndrome()