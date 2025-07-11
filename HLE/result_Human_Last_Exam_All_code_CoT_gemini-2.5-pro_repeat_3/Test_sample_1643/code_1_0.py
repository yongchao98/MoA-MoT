import sys

def analyze_aps_diagnosis():
    """
    Analyzes patient data to determine if they meet the criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    patient_clinical = {
        "vte_events": 3,
    }

    lab_t1 = {  # 3 months ago
        "aCL_IgM": 32,
        "aCL_IgG": 9,
        "b2gp1_IgM": 41,
        "b2gp1_IgG": 18,
        "dRVVT_ratio": 1.44,
    }

    lab_t2 = {  # Today
        "aCL_IgM": 47,
        "aCL_IgG": 7,
        "b2gp1_IgM": 29,
        "b2gp1_IgG": 21,
        "dRVVT_ratio": 1.51,
    }

    # --- APS Criteria Thresholds ---
    criteria = {
        "vte_min_events": 1,
        "aCL_norm": 20,
        "aCL_medium_high_titer": 40,
        "b2gp1_norm": 20,
        "dRVVT_norm": 1.2,
    }

    print("Analyzing patient for Antiphospholipid Syndrome (APS) Diagnosis")
    print("=" * 60)

    # --- Step 1: Check Clinical Criteria ---
    print("Step 1: Evaluating Clinical Criteria")
    clinical_criterion_met = False
    if patient_clinical["vte_events"] >= criteria["vte_min_events"]:
        clinical_criterion_met = True
        print(f"  - Result: MET")
        print(f"  - Reason: Patient has {patient_clinical['vte_events']} VTE events, which is >= the required {criteria['vte_min_events']} event(s).")
    else:
        print("  - Result: NOT MET")
    print("-" * 60)

    # --- Step 2: Check Laboratory Criteria ---
    print("Step 2: Evaluating Laboratory Criteria (Positivity on 2 occasions >= 12 weeks apart)")
    lab_criteria_met_count = 0

    # 2a: Lupus Anticoagulant (LA) via dRVVT
    print("  2a) Lupus Anticoagulant (LA):")
    la_t1_pos = lab_t1["dRVVT_ratio"] > criteria["dRVVT_norm"]
    la_t2_pos = lab_t2["dRVVT_ratio"] > criteria["dRVVT_norm"]
    if la_t1_pos and la_t2_pos:
        lab_criteria_met_count += 1
        print("    - Result: MET")
        print(f"    - Reason: Persistently positive dRVVT ratio.")
        print(f"      - Test 1 dRVVT Ratio: {lab_t1['dRVVT_ratio']} (Normal < {criteria['dRVVT_norm']})")
        print(f"      - Test 2 dRVVT Ratio: {lab_t2['dRVVT_ratio']} (Normal < {criteria['dRVVT_norm']})")
    else:
        print("    - Result: NOT MET")

    # 2b: Anticardiolipin (aCL) Antibody
    print("\n  2b) Anticardiolipin (aCL) Antibody:")
    acl_t1_pos = lab_t1["aCL_IgM"] > criteria["aCL_norm"] or lab_t1["aCL_IgG"] > criteria["aCL_norm"]
    acl_t2_pos = lab_t2["aCL_IgM"] > criteria["aCL_norm"] or lab_t2["aCL_IgG"] > criteria["aCL_norm"]
    acl_high_titer = (
        lab_t1["aCL_IgM"] > criteria["aCL_medium_high_titer"] or
        lab_t1["aCL_IgG"] > criteria["aCL_medium_high_titer"] or
        lab_t2["aCL_IgM"] > criteria["aCL_medium_high_titer"] or
        lab_t2["aCL_IgG"] > criteria["aCL_medium_high_titer"]
    )
    if acl_t1_pos and acl_t2_pos and acl_high_titer:
        lab_criteria_met_count += 1
        print("    - Result: MET")
        print("    - Reason: Persistently positive aCL antibody with medium-high titer.")
        print(f"      - Test 1 aCL IgM: {lab_t1['aCL_IgM']} UI/L (Normal < {criteria['aCL_norm']})")
        print(f"      - Test 2 aCL IgM: {lab_t2['aCL_IgM']} UI/L (Normal < {criteria['aCL_norm']}; Titer > {criteria['aCL_medium_high_titer']})")
    else:
        print("    - Result: NOT MET")

    # 2c: Anti-β2-glycoprotein-I (anti-β2GPI) Antibody
    print("\n  2c) Anti-β2-glycoprotein-I Antibody:")
    b2gp1_t1_pos = lab_t1["b2gp1_IgM"] > criteria["b2gp1_norm"] or lab_t1["b2gp1_IgG"] > criteria["b2gp1_norm"]
    b2gp1_t2_pos = lab_t2["b2gp1_IgM"] > criteria["b2gp1_norm"] or lab_t2["b2gp1_IgG"] > criteria["b2gp1_norm"]
    if b2gp1_t1_pos and b2gp1_t2_pos:
        lab_criteria_met_count += 1
        print("    - Result: MET")
        print("    - Reason: Persistently positive Anti-β2GP1 antibody.")
        print(f"      - Test 1 Anti-β2GP1 IgM: {lab_t1['b2gp1_IgM']} UI/L (Normal < {criteria['b2gp1_norm']})")
        print(f"      - Test 2 Anti-β2GP1 IgM: {lab_t2['b2gp1_IgM']} UI/L (Normal < {criteria['b2gp1_norm']})")
        print(f"      - Test 2 Anti-β2GP1 IgG: {lab_t2['b2gp1_IgG']} UI/L (Normal < {criteria['b2gp1_norm']})")
    else:
        print("    - Result: NOT MET")

    print("-" * 60)

    # --- Step 3: Final Conclusion ---
    print("Step 3: Final Conclusion")
    if clinical_criterion_met and lab_criteria_met_count > 0:
        final_answer = "Yes"
        print("The patient meets at least one clinical AND at least one laboratory criterion.")
        print("This patient categorizes as having Antiphospholipid Syndrome.")
    else:
        final_answer = "No"
        print("The patient does not meet the criteria for a definite diagnosis of Antiphospholipid Syndrome.")

    # Using sys.stdout to ensure the final answer is on a new line after the main text
    sys.stdout.write("\n")
    sys.stdout.flush()
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    analyze_aps_diagnosis()