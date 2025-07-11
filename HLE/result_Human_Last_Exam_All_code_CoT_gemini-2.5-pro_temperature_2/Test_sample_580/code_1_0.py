def solve_clinical_case():
    """
    This function analyzes the clinical vignette to determine the correct
    diagnostic maneuver.
    """

    # Clinical Analysis
    patient_symptoms = "Pain in the lower right extremity L4-S1 distribution (sciatica)."
    imaging_results = "X-ray imaging is unremarkable."
    likely_diagnosis = "Piriformis Syndrome, where the piriformis muscle compresses the sciatic nerve."

    # Rationale for the diagnostic maneuver
    piriformis_muscle_function = "The primary function of the piriformis muscle is external rotation of the hip."
    provocative_test_logic = (
        "To confirm Piriformis Syndrome, a physician will perform a test that "
        "causes the piriformis muscle to contract against resistance. "
        "This action will compress the sciatic nerve and reproduce the patient's pain."
    )

    # Evaluating the options
    correct_action = "External Rotation"
    correct_option = "D"

    # Print the step-by-step reasoning
    print("Patient analysis suggests a diagnosis of Piriformis Syndrome due to sciatica-like symptoms with a normal X-ray.")
    print("-" * 20)
    print(f"The function of the key muscle involved (Piriformis) is: {piriformis_muscle_function}.")
    print(f"The goal of the test is to have the patient contract this muscle against resistance.")
    print(f"Therefore, the physician will resist the action of: {correct_action}.")
    print("-" * 20)
    print(f"The action that will confirm the diagnosis is: {correct_option}. {correct_action}")

solve_clinical_case()