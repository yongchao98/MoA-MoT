import sys

def analyze_aps_criteria():
    """
    Analyzes patient data to determine if they meet the criteria for Antiphospholipid Syndrome (APS).
    The diagnosis is based on the 2006 revised Sydney classification criteria.
    """

    # --- Patient Data ---
    # Clinical Information
    vte_events = 3
    # Laboratory Results (Timespan between tests is ~3 months, meeting the >= 12 weeks requirement)
    lab_set_1 = {
        'dRVVT_ratio': 1.44,       # N < 1.2
        'anticardiolipin_IgM': 32, # N < 20
        'anti_b2gp1_IgM': 41       # N < 20
    }
    lab_set_2 = {
        'dRVVT_ratio': 1.51,       # N < 1.2
        'anticardiolipin_IgM': 47, # N < 20
        'anti_b2gp1_IgM': 29       # N < 20
    }

    # --- Diagnostic Thresholds ---
    threshold = {
        'dRVVT_ratio': 1.2,
        'antibody_positive': 20,
        'antibody_medium_high_titer': 40
    }

    # --- Analysis ---
    print("Evaluating patient for Antiphospholipid Syndrome (APS) using Sydney classification criteria.\n")

    # 1. Assess Clinical Criteria
    print("Step 1: Assessing Clinical Criteria...")
    clinical_criterion_met = False
    if vte_events >= 1:
        print(f"-> Clinical criterion MET. Patient has a history of {vte_events} vascular thrombosis events.")
        clinical_criterion_met = True
    else:
        print("-> Clinical criterion NOT MET.")

    print("-" * 30)

    # 2. Assess Laboratory Criteria (must be persistent over >= 12 weeks)
    print("Step 2: Assessing Laboratory Criteria...")
    laboratory_criteria_met = False
    met_labs = []

    # Check for persistent Lupus Anticoagulant (via dRVVT)
    if lab_set_1['dRVVT_ratio'] > threshold['dRVVT_ratio'] and lab_set_2['dRVVT_ratio'] > threshold['dRVVT_ratio']:
        laboratory_criteria_met = True
        met_labs.append("Lupus Anticoagulant (LA)")
        print(f"-> Persistent LA (dRVVT) MET:")
        print(f"   Test 1: {lab_set_1['dRVVT_ratio']} (Normal < {threshold['dRVVT_ratio']})")
        print(f"   Test 2: {lab_set_2['dRVVT_ratio']} (Normal < {threshold['dRVVT_ratio']})")

    # Check for persistent anticardiolipin antibodies
    if lab_set_1['anticardiolipin_IgM'] > threshold['antibody_positive'] and lab_set_2['anticardiolipin_IgM'] > threshold['antibody_positive']:
        # Criteria also require medium-high titer (>40 or >99th percentile)
        if lab_set_1['anticardiolipin_IgM'] > threshold['antibody_medium_high_titer'] or lab_set_2['anticardiolipin_IgM'] > threshold['antibody_medium_high_titer']:
            laboratory_criteria_met = True
            if "Anticardiolipin Antibody" not in met_labs:
                met_labs.append("Anticardiolipin Antibody")
            print(f"-> Persistent anticardiolipin IgM MET:")
            print(f"   Test 1: {lab_set_1['anticardiolipin_IgM']} UI/L (Normal < {threshold['antibody_positive']})")
            print(f"   Test 2: {lab_set_2['anticardiolipin_IgM']} UI/L (Normal < {threshold['antibody_positive']}) -> Medium titer")

    # Check for persistent anti-B2GP1 antibodies
    if lab_set_1['anti_b2gp1_IgM'] > threshold['antibody_positive'] and lab_set_2['anti_b2gp1_IgM'] > threshold['antibody_positive']:
        laboratory_criteria_met = True
        if "Anti-B2GP1 Antibody" not in met_labs:
            met_labs.append("Anti-B2GP1 Antibody")
        print(f"-> Persistent anti-B2GP1 IgM MET:")
        print(f"   Test 1: {lab_set_1['anti_b2gp1_IgM']} UI/L (Normal < {threshold['antibody_positive']})")
        print(f"   Test 2: {lab_set_2['anti_b2gp1_IgM']} UI/L (Normal < {threshold['antibody_positive']})")

    if not laboratory_criteria_met:
        print("-> Laboratory criteria NOT MET.")
    else:
        print(f"\nLaboratory criteria MET based on persistent positive results for: {', '.join(met_labs)}.")

    print("-" * 30)

    # 3. Final Conclusion
    print("Step 3: Final Diagnosis...")
    final_answer = ""
    if clinical_criterion_met and laboratory_criteria_met:
        final_answer = "Yes"
        print(f"Conclusion: The patient meets both clinical and laboratory criteria.")
    else:
        final_answer = "No"
        print("Conclusion: The patient does not meet both clinical and laboratory criteria.")
    
    # Required output format
    sys.stdout.flush()
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    analyze_aps_criteria()