def diagnose_leg_pain(patient_info, exam_findings, test_options):
    """
    Analyzes clinical information to determine the most likely diagnosis and confirmatory test.

    Args:
        patient_info (dict): A dictionary of the patient's key data.
        exam_findings (list): A list of key examination findings.
        test_options (dict): A dictionary of possible physical exam maneuvers.

    Returns:
        str: A formatted string explaining the diagnosis and the correct action.
    """
    # Step 1: Analyze symptoms and rule out common spinal causes.
    # Sciatica symptoms with an unremarkable X-ray points towards a soft-tissue cause.
    if "sciatica symptoms" in exam_findings and "unremarkable x-ray" in exam_findings:
        likely_diagnosis = "Piriformis Syndrome"
    else:
        likely_diagnosis = "Undetermined"

    # Step 2: Determine the confirmatory test based on the likely diagnosis.
    # Piriformis syndrome is caused by the piriformis muscle compressing the sciatic nerve.
    # The piriformis muscle is a primary EXTERNAL rotator of the hip.
    # To confirm the diagnosis, the muscle must be stretched to reproduce symptoms.
    # Stretching an external rotator is achieved by performing the opposite motion: INTERNAL rotation.
    piriformis_function = "External Rotation"
    confirmatory_action = "Internal Rotation"

    # Step 3: Find the correct answer choice.
    correct_choice = None
    for key, value in test_options.items():
        if value == confirmatory_action:
            correct_choice = key
            break

    # Step 4: Formulate the final answer, including key numbers from the case.
    # The prompt requests including numbers in the final output.
    output = (
        f"Patient Profile: Age {patient_info['age']}, BMI {patient_info['bmi']}, Pain History {patient_info['pain_history_months']} months.\n"
        f"Diagnosis Logic:\n"
        f"1. The patient's sciatica-like symptoms with an unremarkable X-ray strongly suggest a soft-tissue cause like {likely_diagnosis}.\n"
        f"2. The piriformis muscle's primary function is {piriformis_function} of the hip.\n"
        f"3. To provoke symptoms and confirm the diagnosis, the muscle must be stretched by performing the opposing action.\n"
        f"4. The opposite of {piriformis_function} is {confirmatory_action}.\n"
        f"Conclusion: The action that will confirm the diagnosis is {confirmatory_action}.\n"
        f"Final Answer Choice: {correct_choice}"
    )

    return output

# Patient data from the clinical case
patient_data = {
    "age": 55,
    "bmi": 29,
    "pain_history_months": 5
}

# Key findings from the exam
clinical_findings = [
    "sciatica symptoms", # Pain in L4-S1 distribution
    "unremarkable x-ray",
    "pain intensified by lying supine",
    "test position: left decubitus"
]

# Provided answer choices
maneuver_options = {
    "A": "Abduction",
    "B": "Adduction",
    "C": "Internal Rotation",
    "D": "External Rotation",
    "E": "Flexion",
    "F": "Extension"
}

# Run the diagnostic logic and print the result
final_explanation = diagnose_leg_pain(patient_data, clinical_findings, maneuver_options)
print(final_explanation)