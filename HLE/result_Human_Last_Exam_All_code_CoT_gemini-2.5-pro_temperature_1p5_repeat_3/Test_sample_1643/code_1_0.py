import sys

# Suppress warnings for cleaner output, not essential for logic
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def diagnose_aps():
    """
    Analyzes patient data to determine if they meet the criteria for
    Antiphospholipid Syndrome (APS) based on the revised Sapporo criteria.
    """

    # --- Patient's Clinical Data ---
    # The primary clinical criterion is vascular thrombosis.
    vte_event_count = 3  # Calf DVT, PE, another PE

    # --- Patient's Laboratory Data ---
    # Measured at two time points, 3 months (~12 weeks) apart.
    # Timepoint 1 (3 months ago)
    lab1_anticardiolipin_IgM = 32  # Normal < 20 UI/L
    lab1_antiB2GP1_IgM = 41         # Normal < 20 UI/L
    lab1_dRVVT_ratio = 1.44          # Normal < 1.2

    # Timepoint 2 (Today)
    lab2_anticardiolipin_IgM = 47  # Normal < 20 UI/L
    lab2_antiB2GP1_IgM = 29         # Normal < 20 UI/L
    lab2_dRVVT_ratio = 1.51          # Normal < 1.2

    # --- Diagnostic Logic ---
    clinical_criterion_met = False
    lab_criteria_met = {
        'Lupus Anticoagulant': False,
        'Anticardiolipin': False,
        'AntiB2GP1': False
    }

    # Step 1: Evaluate Clinical Criteria
    print("--- Evaluating Clinical Criteria ---")
    if vte_event_count >= 1:
        clinical_criterion_met = True
        print(f"Result: MET. Patient has a history of {vte_event_count} venous thrombotic events.")
    else:
        print(f"Result: NOT MET. Patient has no documented thrombosis.")

    # Step 2: Evaluate Laboratory Criteria (must be persistent over 12 weeks)
    print("\n--- Evaluating Laboratory Criteria (Persistence over >12 weeks) ---")

    # A) Lupus Anticoagulant (LA), assessed by dRVVT
    is_la_positive_t1 = lab1_dRVVT_ratio > 1.2
    is_la_positive_t2 = lab2_dRVVT_ratio > 1.2
    if is_la_positive_t1 and is_la_positive_t2:
        lab_criteria_met['Lupus Anticoagulant'] = True
        print(f"1. Lupus Anticoagulant: MET")
        print(f"   - Positive 3 months ago (dRVVT Ratio: {lab1_dRVVT_ratio}) and today (dRVVT Ratio: {lab2_dRVVT_ratio}).")
    else:
        print("1. Lupus Anticoagulant: NOT MET")

    # B) Anticardiolipin (aCL) Antibodies (IgM or IgG)
    # Criterion requires medium-to-high titer (>40 units, or >99th percentile).
    is_acl_positive_t1 = lab1_anticardiolipin_IgM > 20
    is_acl_positive_t2 = lab2_anticardiolipin_IgM > 20
    has_medium_high_titer_acl = lab1_anticardiolipin_IgM > 40 or lab2_anticardiolipin_IgM > 40
    if is_acl_positive_t1 and is_acl_positive_t2 and has_medium_high_titer_acl:
        lab_criteria_met['Anticardiolipin'] = True
        print(f"2. Anticardiolipin IgM Antibody: MET")
        print(f"   - Persistently positive with a value of {lab2_anticardiolipin_IgM} meeting medium-titer criteria (>40).")
    else:
        print("2. Anticardiolipin IgM Antibody: NOT MET")

    # C) Anti-β2-Glycoprotein-I (aβ2GPI) Antibodies (IgM or IgG)
    is_ab2gp1_positive_t1 = lab1_antiB2GP1_IgM > 20
    is_ab2gp1_positive_t2 = lab2_antiB2GP1_IgM > 20
    if is_ab2gp1_positive_t1 and is_ab2gp1_positive_t2:
        lab_criteria_met['AntiB2GP1'] = True
        print(f"3. Anti-β2GP1 IgM Antibody: MET")
        print(f"   - Persistently positive (Values: {lab1_antiB2GP1_IgM} and {lab2_antiB2GP1_IgM}).")
    else:
        print("3. Anti-β2GP1 IgM Antibody: NOT MET")
    
    # Check if at least one laboratory criterion is fulfilled
    any_lab_criterion_met = any(lab_criteria_met.values())

    # Step 3: Final Conclusion
    print("\n--- Final Diagnosis ---")
    if clinical_criterion_met and any_lab_criterion_met:
        final_answer = "Yes"
        print("CONCLUSION: Patient meets at least one clinical criterion AND at least one laboratory criterion.")
    else:
        final_answer = "No"
        print("CONCLUSION: Patient does not meet the full diagnostic criteria for APS.")

    # Final formatted answer
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"\n<<<{final_answer}>>>")

if __name__ == "__main__":
    diagnose_aps()