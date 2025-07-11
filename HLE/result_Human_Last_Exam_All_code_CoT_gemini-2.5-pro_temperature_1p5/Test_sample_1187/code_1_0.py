def find_lesion_level():
    """
    Identifies the spinal cord lesion level based on anatomical landmarks.
    """

    # A dictionary mapping key dermatome levels to their anatomical landmarks.
    dermatome_map = {
        'T4': 'Nipple line',
        'T6': 'Xiphoid process',
        'T10': 'Umbilicus (belly button)',
        'T12': 'Inguinal region',
        'L4': 'Knee'
    }

    # The patient's sensory deficit begins at the umbilicus.
    sensory_level_landmark = 'Umbilicus (belly button)'

    # Find the corresponding spinal cord level.
    injury_level = None
    for level, landmark in dermatome_map.items():
        if landmark == sensory_level_landmark:
            injury_level = level
            break

    print("Step 1: The patient's symptoms indicate a Brown-Séquard syndrome (spinal cord hemisection).")
    print("Step 2: The level of injury is determined by the sensory deficit.")
    print(f"Step 3: The sensory loss for pain and temperature starts at the '{sensory_level_landmark}'.")
    print(f"Step 4: The anatomical landmark '{sensory_level_landmark}' corresponds to the {injury_level} dermatome.")
    print(f"Step 5: Therefore, the spinal cord injury is located at the {injury_level} level.")

find_lesion_level()