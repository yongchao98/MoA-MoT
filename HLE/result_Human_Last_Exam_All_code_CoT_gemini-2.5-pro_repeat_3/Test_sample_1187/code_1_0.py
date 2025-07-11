def find_injury_location():
    """
    This function analyzes the clinical presentation to determine the spinal cord injury level.
    """
    # Key clinical finding for localization is the sensory level.
    sensory_loss_landmark = "umbilicus"
    
    # Neurological mapping of dermatomes to spinal levels.
    dermatome_map = {
        "C7": "Middle finger",
        "T4": "Nipple line",
        "T6": "Xiphoid process",
        "T10": "Umbilicus (belly button)",
        "T12": "Inguinal or groin region",
        "L4": "Knee",
        "L5": "Top of the foot and big toe"
    }

    # The patient presents with classic signs of Brown-Séquard syndrome.
    # The key to locating the lesion vertically is the sensory level for pain and temperature.
    print("Patient Symptom Analysis:")
    print("1. Ipsilateral (right side) motor weakness and loss of vibration/proprioception.")
    print("2. Contralateral (left side) loss of pain and temperature sensation.")
    print("This pattern indicates Brown-Séquard Syndrome (a spinal cord hemisection).\n")

    print("Determining the Spinal Level:")
    print(f"The sensory loss for pain and temperature is described as starting at the '{sensory_loss_landmark}'.")

    # Find the corresponding spinal level for the landmark.
    injury_level = None
    for level, landmark in dermatome_map.items():
        if sensory_loss_landmark in landmark:
            injury_level = level
            break
            
    if injury_level:
        print(f"The dermatome corresponding to the {sensory_loss_landmark} is {injury_level}.")
        print(f"Therefore, the spinal cord injury is located at the {injury_level} level.")
    else:
        print("Could not determine the injury level from the provided landmark.")

find_injury_location()
<<<H>>>