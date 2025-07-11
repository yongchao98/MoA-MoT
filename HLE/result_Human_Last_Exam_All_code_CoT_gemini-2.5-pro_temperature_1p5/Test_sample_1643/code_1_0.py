def check_antiphospholipid_syndrome():
    """
    Analyzes patient data against the Sydney criteria for Antiphospholipid Syndrome (APS).
    A diagnosis requires at least one clinical and one laboratory criterion.
    """

    # --- Patient Data ---
    vte_events = 3
    # Lab data from 3 months ago
    lab1 = {
        "dRVVT_ratio": 1.44, "dRVVT_norm": 1.2,
        "aCL_IgM": 32, "aCL_norm": 20,
        "aCL_IgG": 9,
        "antiB2GP1_IgM": 41, "antiB2GP1_norm": 20,
        "antiB2GP1_IgG": 18,
    }
    # Lab data from today (interval > 12 weeks)
    lab2 = {
        "dRVVT_ratio": 1.51, "dRVVT_norm": 1.2,
        "aCL_IgM": 47, "aCL_norm": 20,
        "aCL_IgG": 7,
        "antiB2GP1_IgM": 29, "antiB2GP1_norm": 20,
        "antiB2GP1_IgG": 21,
    }

    # --- Step 1: Evaluate Clinical Criteria ---
    print("--- Evaluating Clinical Criteria for APS ---")
    has_clinical_criterion = vte_events >= 1
    if has_clinical_criterion:
        print(f"Clinical Criterion Met: Yes (Patient has {vte_events} documented VTE events).")
    else:
        print("Clinical Criterion Met: No.")

    # --- Step 2: Evaluate Laboratory Criteria (Persistence over >12 weeks) ---
    print("\n--- Evaluating Laboratory Criteria for APS ---")

    # Check for persistent Lupus Anticoagulant (LA)
    persistent_la = lab1["dRVVT_ratio"] > lab1["dRVVT_norm"] and lab2["dRVVT_ratio"] > lab2["dRVVT_norm"]
    print(f"1. Persistent Lupus Anticoagulant (dRVVT > {lab1['dRVVT_norm']})?")
    print(f"   - 3 months ago: {lab1['dRVVT_ratio']} | Today: {lab2['dRVVT_ratio']} | Met: {persistent_la}")

    # Check for persistent Anticardiolipin (aCL) antibodies
    persistent_acl = (lab1["aCL_IgM"] > lab1["aCL_norm"] and lab2["aCL_IgM"] > lab2["aCL_norm"]) or \
                     (lab1["aCL_IgG"] > lab1["aCL_norm"] and lab2["aCL_IgG"] > lab2["aCL_norm"])
    print(f"2. Persistent Anticardiolipin Ab (aCL > {lab1['aCL_norm']})?")
    print(f"   - IgM 3mo: {lab1['aCL_IgM']} | IgM Today: {lab2['aCL_IgM']} | IgG 3mo: {lab1['aCL_IgG']} | IgG Today: {lab2['aCL_IgG']} | Met: {persistent_acl}")


    # Check for persistent Anti-β2 Glycoprotein I (anti-B2GP1) antibodies
    persistent_b2gp1 = (lab1["antiB2GP1_IgM"] > lab1["antiB2GP1_norm"] and lab2["antiB2GP1_IgM"] > lab2["antiB2GP1_norm"]) or \
                       (lab1["antiB2GP1_IgG"] > lab1["antiB2GP1_norm"] and lab2["antiB2GP1_IgG"] > lab2["antiB2GP1_norm"])
    print(f"3. Persistent Anti-B2GP1 Ab (anti-B2GP1 > {lab1['antiB2GP1_norm']})?")
    print(f"   - IgM 3mo: {lab1['antiB2GP1_IgM']} | IgM Today: {lab2['antiB2GP1_IgM']} | IgG 3mo: {lab1['antiB2GP1_IgG']} | IgG Today: {lab2['antiB2GP1_IgG']} | Met: {persistent_b2gp1}")


    has_lab_criterion = persistent_la or persistent_acl or persistent_b2gp1
    if has_lab_criterion:
        print("\nAt least one laboratory criterion has been met.")
    else:
        print("\nNo laboratory criteria have been met.")

    # --- Step 3: Final Conclusion ---
    print("\n--- Final Diagnosis ---")
    print(f"Clinical Criterion Met: {has_clinical_criterion}")
    print(f"Laboratory Criterion Met: {has_lab_criterion}")
    
    if has_clinical_criterion and has_lab_criterion:
        final_answer = "Yes"
    else:
        final_answer = "No"

    print("\nBased on the analysis, does this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")

if __name__ == "__main__":
    check_antiphospholipid_syndrome()