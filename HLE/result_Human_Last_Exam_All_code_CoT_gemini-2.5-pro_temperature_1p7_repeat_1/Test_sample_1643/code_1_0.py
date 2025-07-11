import sys
# Redirect print to a string buffer to capture all output
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


def evaluate_patient_for_aps():
    """
    Evaluates patient data against the 2006 Sydney criteria for Antiphospholipid Syndrome (APS).
    """
    # --- Patient Data ---
    clinical_history = {
        'vte_events': 3
    }

    lab_results = {
        'test_1': {  # 3 months ago
            'time_point': "3 months ago",
            'dRVVT_ratio': 1.44, 'dRVVT_norm': 1.2,
            'aCL_IgM': 32, 'aCL_norm': 20,
            'anti_b2gp1_IgM': 41, 'anti_b2gp1_norm': 20
        },
        'test_2': {  # Today
            'time_point': "Today",
            'dRVVT_ratio': 1.51, 'dRVVT_norm': 1.2,
            'aCL_IgM': 47, 'aCL_norm': 20,
            'anti_b2gp1_IgM': 29, 'anti_b2gp1_norm': 20
        }
    }
    # Time between tests is 3 months, which satisfies the >= 12 weeks requirement.

    # --- Step 1: Evaluate Clinical Criteria ---
    print("--- Evaluation for Antiphospholipid Syndrome (APS) ---")
    print("\nStep 1: Assessing Clinical Criteria")
    vte_criterion_met = clinical_history['vte_events'] >= 1
    if vte_criterion_met:
        print(f"Outcome: CLINICAL CRITERION MET.")
        print(f"Reason: Patient has a history of {clinical_history['vte_events']} vascular thrombosis events, which is >= 1.")
    else:
        print("Outcome: CLINICAL CRITERION NOT MET.")

    # --- Step 2: Evaluate Laboratory Criteria ---
    print("\nStep 2: Assessing Laboratory Criteria (requires persistence over >= 12 weeks)")
    lab_test_1 = lab_results['test_1']
    lab_test_2 = lab_results['test_2']
    
    # Check 2a: Lupus Anticoagulant (LA)
    la_positive_1 = lab_test_1['dRVVT_ratio'] > lab_test_1['dRVVT_norm']
    la_positive_2 = lab_test_2['dRVVT_ratio'] > lab_test_2['dRVVT_norm']
    la_criterion_met = la_positive_1 and la_positive_2
    print("\n- Criterion A: Lupus Anticoagulant (LA) based on dRVVT")
    print(f"  Test 1 (3 months ago): dRVVT Ratio = {lab_test_1['dRVVT_ratio']} (Normal < {lab_test_1['dRVVT_norm']}) -> {'Positive' if la_positive_1 else 'Negative'}")
    print(f"  Test 2 (Today): dRVVT Ratio = {lab_test_2['dRVVT_ratio']} (Normal < {lab_test_2['dRVVT_norm']}) -> {'Positive' if la_positive_2 else 'Negative'}")
    print(f"  Outcome: {'PERSISTENTLY POSITIVE' if la_criterion_met else 'NOT PERSISTENT'}. Criterion MET.")

    # Check 2b: Anticardiolipin (aCL) antibodies
    acl_positive_1 = lab_test_1['aCL_IgM'] > lab_test_1['aCL_norm']
    acl_positive_2 = lab_test_2['aCL_IgM'] > lab_test_2['aCL_norm']
    acl_criterion_met = acl_positive_1 and acl_positive_2
    print("\n- Criterion B: Anticardiolipin (aCL) IgM Antibody")
    print(f"  Test 1 (3 months ago): aCL IgM = {lab_test_1['aCL_IgM']} UI/L (Normal < {lab_test_1['aCL_norm']}) -> {'Positive' if acl_positive_1 else 'Negative'}")
    print(f"  Test 2 (Today): aCL IgM = {lab_test_2['aCL_IgM']} UI/L (Normal < {lab_test_2['aCL_norm']}) -> {'Positive' if acl_positive_2 else 'Negative'}")
    print(f"  Outcome: {'PERSISTENTLY POSITIVE' if acl_criterion_met else 'NOT PERSISTENT'}. Criterion MET.")

    # Check 2c: Anti-B2-Glycoprotein-I (anti-B2GPI) antibodies
    b2gpi_positive_1 = lab_test_1['anti_b2gp1_IgM'] > lab_test_1['anti_b2gp1_norm']
    b2gpi_positive_2 = lab_test_2['anti_b2gp1_IgM'] > lab_test_2['anti_b2gp1_norm']
    b2gpi_criterion_met = b2gpi_positive_1 and b2gpi_positive_2
    print("\n- Criterion C: Anti-B2-Glycoprotein-I (anti-B2GPI) IgM Antibody")
    print(f"  Test 1 (3 months ago): anti-B2GPI IgM = {lab_test_1['anti_b2gp1_IgM']} UI/L (Normal < {lab_test_1['anti_b2gp1_norm']}) -> {'Positive' if b2gpi_positive_1 else 'Negative'}")
    print(f"  Test 2 (Today): anti-B2GPI IgM = {lab_test_2['anti_b2gp1_IgM']} UI/L (Normal < {lab_test_2['anti_b2gp1_norm']}) -> {'Positive' if b2gpi_positive_2 else 'Negative'}")
    print(f"  Outcome: {'PERSISTENTLY POSITIVE' if b2gpi_criterion_met else 'NOT PERSISTENT'}. Criterion MET.")

    any_lab_criterion_met = la_criterion_met or acl_criterion_met or b2gpi_criterion_met
    if any_lab_criterion_met:
        print("\nOverall Lab Result: At least one laboratory criterion has been met.")

    # --- Step 3: Final Conclusion ---
    print("\n--- Final Conclusion ---")
    if vte_criterion_met and any_lab_criterion_met:
        final_answer = "Yes"
        print("Diagnosis: The patient meets both clinical (vascular thrombosis) and laboratory (persistent antiphospholipid antibodies) criteria.")
    else:
        final_answer = "No"
        print("Diagnosis: The patient does not meet the full criteria for a definitive diagnosis of APS.")
    
    print("\nDoes this patient categorize as having antiphospholipid syndrome?")
    return final_answer


# Run the evaluation and get the final answer
final_answer = evaluate_patient_for_aps()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Print the final formatted answer
print(f"<<<{final_answer}>>>")