def find_injury_location():
    """
    Analyzes clinical findings to determine the location of a spinal cord injury.
    """
    # Patient's key neurological findings
    ipsilateral_motor_loss = "Right leg weakness"
    ipsilateral_proprioception_vibration_loss = "Right leg"
    contralateral_pain_temp_loss = "Left side"
    sensory_level_landmark = "Umbilicus"

    # Dermatome mapping
    dermatome_map = {
        "Nipple line": "T4",
        "Xiphoid process": "T6",
        "Umbilicus": "T10",
        "Inguinal region": "T12"
    }

    print("Analyzing the neurological presentation:")
    print(f"- Ipsilateral (Right Side) Signs: {ipsilateral_motor_loss} and loss of proprioception/vibration.")
    print(f"- Contralateral (Left Side) Signs: Loss of pain and temperature sensation.")
    print("\nThis pattern is characteristic of Brown-Séquard Syndrome (a spinal cord hemisection).")

    # Determine the level of injury from the sensory landmark
    lesion_level_dermatome = dermatome_map.get(sensory_level_landmark)

    print(f"\nThe sensory level for pain and temperature loss is at the {sensory_level_landmark}.")
    print(f"The dermatome corresponding to the {sensory_level_landmark} is {lesion_level_dermatome}.")
    
    # Final conclusion output
    print("\n--- Final Diagnosis ---")
    print(f"The lesion is a hemisection of the spinal cord located at the {lesion_level_dermatome} level.")
    print("Therefore, the corresponding answer choice is H.")


find_injury_location()
<<<H>>>