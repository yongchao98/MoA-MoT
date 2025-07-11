def diagnose_patient_condition():
    """
    This program models the clinical reasoning used to answer the user's question.
    It analyzes the patient's symptoms and the described physical exam maneuver
    to determine the most likely diagnosis and the specific action that would confirm it.
    """

    # Patient details and symptoms suggest sciatic nerve involvement without bony pathology.
    leading_differential = "Piriformis Syndrome"
    
    # The physical exam is performed with the patient in the left decubitus position,
    # testing the extended right leg against resistance.
    exam_setup = {
        "position": "Left decubitus (lying on left side)",
        "leg_tested": "Right",
        "leg_position": "Extended",
        "action": "Applying resistance"
    }

    # The piriformis muscle is a primary external rotator of the hip when it is extended.
    # In Piriformis Syndrome, the muscle compresses the sciatic nerve.
    # Therefore, making the muscle contract against resistance should reproduce the pain.
    
    actions = {
        "A": "Abduction",
        "B": "Adduction",
        "C": "Internal Rotation",
        "D": "External Rotation",
        "E": "Flexion",
        "F": "Extension"
    }
    
    reasoning = {
        "A": "Tests gluteus medius/minimus, not the primary muscle of interest.",
        "B": "Tests adductor muscles.",
        "C": "Stretches the piriformis, but is not a test of its contraction.",
        "D": "Directly causes contraction of the piriformis muscle against resistance, a classic test for Piriformis Syndrome.",
        "E": "Tests hip flexors.",
        "F": "Tests gluteus maximus."
    }

    print("Analyzing the Clinical Scenario:")
    print(f"1. The patient's symptoms point towards a sciatic nerve issue, such as {leading_differential}.")
    print(f"2. The physical exam setup ({exam_setup['position']}, testing the {exam_setup['leg_position']} {exam_setup['leg_tested']} leg) is ideal for testing hip rotator muscles.")
    print("3. The goal is to perform an action that reproduces the pain by stressing the structure causing it.")
    print("4. The piriformis muscle's primary function in this position is external rotation.")
    print("-" * 20)
    print("Evaluating Answer Choices:")

    correct_choice = None
    for key, action in actions.items():
        if action == "External Rotation":
            conclusion = f"This action directly contracts the piriformis muscle. If this reproduces the patient's sciatic pain, it helps confirm {leading_differential}."
            correct_choice = key
        else:
            conclusion = reasoning[key]
        print(f"Choice {key}. {action}: {conclusion}")

    print("-" * 20)
    print(f"Final Conclusion: The most specific action to confirm the diagnosis in this context is resisted {actions[correct_choice]}.")
    print(f"The correct answer is therefore '{correct_choice}'.")

diagnose_patient_condition()