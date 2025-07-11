def diagnose_spinal_injury_level():
    """
    Analyzes a clinical vignette to determine the level of a spinal cord injury.
    """
    # Step 1: Define key clinical findings from the case.
    ipsilateral_motor_loss_side = "right"
    ipsilateral_proprioception_loss_side = "right"
    contralateral_pain_temp_loss_side = "left"
    sensory_level_landmark = "umbilicus"

    # Step 2: Define a mapping of common dermatome landmarks to spinal levels.
    dermatome_landmarks = {
        "nipple": "T4",
        "xiphoid_process": "T6",
        "umbilicus": "T10",
        "inguinal_ligament": "T12"
    }

    # Step 3: Identify the syndrome based on the pattern of deficits.
    # The pattern is classic Brown-Séquard syndrome (spinal cord hemisection).
    # The injury is on the same side as motor/proprioception loss.
    injury_side = ipsilateral_motor_loss_side

    # Step 4: Determine the spinal level from the sensory landmark.
    # The sensory level for pain/temperature loss is the most reliable indicator.
    if sensory_level_landmark in dermatome_landmarks:
        injury_level = dermatome_landmarks[sensory_level_landmark]
        thoracic_level_number = 10
    else:
        injury_level = "Unknown"
        thoracic_level_number = None

    # Step 5: Print the reasoning and the final answer.
    print("Clinical Reasoning:")
    print(f"1. The patient presents with ipsilateral (right-sided) motor weakness and loss of proprioception/vibration.")
    print(f"2. The patient also has contralateral (left-sided) loss of pain and temperature sensation.")
    print(f"3. This specific pattern indicates Brown-Séquard syndrome, a hemisection of the spinal cord, with the injury on the {injury_side} side.")
    print(f"4. The sensory loss begins at the '{sensory_level_landmark}'. This landmark corresponds to a specific dermatome.")
    print(f"5. The dermatome for the umbilicus is T{thoracic_level_number}.")
    print("\nConclusion:")
    print(f"The location of the patient's injury is at the T{thoracic_level_number} spinal level.")

diagnose_spinal_injury_level()