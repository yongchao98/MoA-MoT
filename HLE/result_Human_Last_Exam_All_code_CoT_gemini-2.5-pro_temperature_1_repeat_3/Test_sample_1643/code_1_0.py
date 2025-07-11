import sys
# Redirect print to a string to control output order
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


def check_aps_diagnosis():
    """
    Evaluates if a patient meets the criteria for Antiphospholipid Syndrome (APS)
    based on the revised Sapporo classification criteria.
    """

    # --- Patient Data ---
    # Clinical History
    vte_events_count = 3

    # Lab tests from 3 months ago
    lab1 = {
        "anti_b2gp1_igm": 41,
        "anticardiolipin_igm": 32,
        "drvv_ratio": 1.44
    }

    # Lab tests from today (12 weeks later)
    lab2 = {
        "anti_b2gp1_igm": 29,
        "anticardiolipin_igm": 47,
        "drvv_ratio": 1.51
    }

    # Time between lab tests in weeks
    time_between_labs_weeks = 12

    print("Step 1: Evaluating Clinical Criteria for APS")
    print("---------------------------------------------")
    print("Criteria: At least one clinical episode of venous, arterial, or small vessel thrombosis.")
    
    clinical_criterion_met = vte_events_count >= 1
    
    if clinical_criterion_met:
        print(f"Result: MET. The patient has a history of {vte_events_count} VTE events.")
    else:
        print(f"Result: NOT MET. The patient has {vte_events_count} VTE events.")

    print("\nStep 2: Evaluating Laboratory Criteria for APS")
    print("----------------------------------------------")
    print(f"Criteria: A positive lab test must be confirmed on two occasions at least {time_between_labs_weeks} weeks apart.")
    
    lab_criterion_met = False
    fulfilled_lab_criteria = []

    # a) Lupus Anticoagulant (LA) based on dRVVT
    # Note: Rivaroxaban significantly interferes with dRVVT, often causing false positives.
    # This criterion is therefore unreliable without specific drug neutralization tests.
    la_positive_1 = lab1['drvv_ratio'] > 1.2
    la_positive_2 = lab2['drvv_ratio'] > 1.2
    print("\n- Checking for persistent Lupus Anticoagulant (dRVVT > 1.2)...")
    print(f"  - Test 1 (3 months ago): dRVVT ratio = {lab1['drvv_ratio']}. Positive: {la_positive_1}")
    print(f"  - Test 2 (Today): dRVVT ratio = {lab2['drvv_ratio']}. Positive: {la_positive_2}")
    print("  - Status: Persistently Positive. However, these results are considered UNRELIABLE for diagnosis due to known interference from Rivaroxaban.")

    # b) Anticardiolipin (aCL) IgM antibody (medium/high titer > 40 UI/L)
    acl_positive_1 = lab1['anticardiolipin_igm'] > 20
    acl_positive_2 = lab2['anticardiolipin_igm'] > 20
    medium_high_titer = lab1['anticardiolipin_igm'] > 40 or lab2['anticardiolipin_igm'] > 40
    
    print("\n- Checking for persistent Anticardiolipin IgM (>20 UI/L) with medium/high titer (>40 UI/L)...")
    print(f"  - Test 1 (3 months ago): aCL IgM = {lab1['anticardiolipin_igm']} UI/L. Positive: {acl_positive_1}")
    print(f"  - Test 2 (Today): aCL IgM = {lab2['anticardiolipin_igm']} UI/L. Positive: {acl_positive_2}")
    
    if acl_positive_1 and acl_positive_2 and medium_high_titer:
        lab_criterion_met = True
        fulfilled_lab_criteria.append("Anticardiolipin IgM")
        print("  - Status: MET. The antibody is persistently positive, with one result in the medium titer range (>40 UI/L).")
    else:
        print("  - Status: NOT MET.")

    # c) Anti-β2 glycoprotein-I (anti-β2GPI) IgM antibody (>20 UI/L)
    b2gpi_positive_1 = lab1['anti_b2gp1_igm'] > 20
    b2gpi_positive_2 = lab2['anti_b2gp1_igm'] > 20
    
    print("\n- Checking for persistent Anti-β2GPI IgM (>20 UI/L)...")
    print(f"  - Test 1 (3 months ago): anti-β2GPI IgM = {lab1['anti_b2gp1_igm']} UI/L. Positive: {b2gpi_positive_1}")
    print(f"  - Test 2 (Today): anti-β2GPI IgM = {lab2['anti_b2gp1_igm']} UI/L. Positive: {b2gpi_positive_2}")
    
    if b2gpi_positive_1 and b2gpi_positive_2:
        lab_criterion_met = True
        if "Anti-β2GPI IgM" not in fulfilled_lab_criteria:
            fulfilled_lab_criteria.append("Anti-β2GPI IgM")
        print("  - Status: MET. The antibody is persistently positive.")
    else:
        print("  - Status: NOT MET.")
    
    print("\n--- Final Conclusion ---")
    
    final_answer = "No"
    if clinical_criterion_met and lab_criterion_met:
        final_answer = "Yes"
        print("The patient meets at least one clinical criterion (Vascular Thrombosis) AND at least one laboratory criterion.")
        print(f"Fulfilled laboratory criteria: {', '.join(fulfilled_lab_criteria)}.")
        print("Therefore, the patient categorizes as having Antiphospholipid Syndrome.")
    else:
        print("The patient does not meet both the clinical and laboratory criteria for APS.")

    # Return the final answer for capture
    return final_answer

# Run the function and capture the answer
final_answer_value = check_aps_diagnosis()

# Restore standard output and print the captured analysis
sys.stdout = old_stdout
print(captured_output.getvalue())

# Print the final answer in the required format
print("\nDoes this patient categorizes as having antiphospholipid syndrome?")
print("Answer by Yes or No")
print(f'<<<{final_answer_value}>>>')