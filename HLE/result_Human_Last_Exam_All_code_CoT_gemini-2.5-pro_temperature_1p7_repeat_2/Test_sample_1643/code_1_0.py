import sys

# Suppress writing __pycache__ files
sys.dont_write_bytecode = True

def evaluate_aps_diagnosis():
    """
    Analyzes patient data against the 2006 Sydney criteria for Antiphospholipid Syndrome (APS).
    """

    # --- Patient Data ---
    clinical_history = {
        "vte_events": 3,
        "first_vte": "Calf DVT at age 18 after a 6-hour flight",
        "second_vte": "Intermediate-risk PE at age 25 during pregnancy",
        "third_vte": "Spontaneous PE 4 months ago"
    }

    lab_test_1 = {
        "time": "3 months ago",
        "aCL_IgM": 32,  # N < 20 UI/L
        "anti_b2gp1_IgM": 41, # N < 20 UI/L
        "drvvT_ratio": 1.44 # N < 1.2
    }

    lab_test_2 = {
        "time": "Today",
        "aCL_IgM": 47,  # N < 20 UI/L
        "anti_b2gp1_IgM": 29, # N < 20 UI/L
        "drvvT_ratio": 1.51 # N < 1.2
    }

    on_anticoagulation = True
    anticoagulant_type = "Rivaroxaban (Direct Factor Xa inhibitor)"
    time_between_tests_weeks = 13 # 3 months

    clinical_criterion_met = False
    lab_criterion_met = False

    print("--- Antiphospholipid Syndrome (APS) Diagnostic Evaluation ---")

    # Step 1: Evaluate Clinical Criteria
    print("\nStep 1: Evaluating Clinical Criteria...")
    if clinical_history["vte_events"] >= 1:
        print(f"Finding: Patient has a history of {clinical_history['vte_events']} venous thromboembolic (VTE) events.")
        print("Conclusion: The clinical criterion (Vascular Thrombosis) is MET.")
        clinical_criterion_met = True
    else:
        print("Conclusion: No definitive clinical criterion for vascular thrombosis met.")

    # Step 2: Evaluate Laboratory Criteria
    print("\nStep 2: Evaluating Laboratory Criteria...")
    print(f"Persistence Check: Tests were performed approximately {time_between_tests_weeks} weeks apart, which is >= 12 weeks. The persistence requirement is met.")

    print("\n- Analyzing Lupus Anticoagulant (LA):")
    print(f"  Result 1 (dRVVT): {lab_test_1['drvvT_ratio']} (Normal < 1.2)")
    print(f"  Result 2 (dRVVT): {lab_test_2['drvvT_ratio']} (Normal < 1.2)")
    if on_anticoagulation:
        print(f"  Important Note: Patient is on {anticoagulant_type}, which strongly interferes with clot-based assays like dRVVT, often causing false positives.")
        print("  Conclusion: LA results are UNINTERPRETABLE and cannot be used for diagnosis.")
    
    print("\n- Analyzing Anticardiolipin (aCL) IgM antibodies:")
    print(f"  Result 1: {lab_test_1['aCL_IgM']} UI/L (Normal < 20)")
    print(f"  Result 2: {lab_test_2['aCL_IgM']} UI/L (Normal < 20)")
    # Check for persistence and titer
    is_persistent = lab_test_1['aCL_IgM'] > 20 and lab_test_2['aCL_IgM'] > 20
    is_medium_high_titer = lab_test_1['aCL_IgM'] > 40 or lab_test_2['aCL_IgM'] > 40
    if is_persistent and is_medium_high_titer:
        print(f"  Conclusion: aCL IgM is persistently positive, and the second result ({lab_test_2['aCL_IgM']} UI/L) meets the medium-titer threshold (>40 UI/L). This laboratory criterion is MET.")
        lab_criterion_met = True

    print("\n- Analyzing Anti-ß2 Glycoprotein I (anti-ß2GPI) IgM antibodies:")
    print(f"  Result 1: {lab_test_1['anti_b2gp1_IgM']} UI/L (Normal < 20)")
    print(f"  Result 2: {lab_test_2['anti_b2gp1_IgM']} UI/L (Normal < 20)")
    if lab_test_1['anti_b2gp1_IgM'] > 20 and lab_test_2['anti_b2gp1_IgM'] > 20:
        print("  Conclusion: Anti-ß2GPI IgM is persistently positive. This laboratory criterion is MET.")
        # This also sets the overall lab criterion to True if the aCL part was somehow missed
        lab_criterion_met = True

    # Step 3: Final Conclusion
    print("\n--- Final Diagnostic Conclusion ---")
    print(f"Clinical Criterion Met: {clinical_criterion_met}")
    print(f"Laboratory Criterion Met: {lab_criterion_met}")
    if clinical_criterion_met and lab_criterion_met:
        print("\nDiagnosis: The patient meets both clinical and laboratory criteria for Antiphospholipid Syndrome.")
        print("<<<Yes>>>")
    else:
        print("\nDiagnosis: The patient does not meet the full criteria for Antiphospholipid Syndrome at this time.")
        print("<<<No>>>")

if __name__ == "__main__":
    evaluate_aps_diagnosis()