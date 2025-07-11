def find_injury_location():
    """
    Determines the spinal cord injury level by mapping the key clinical
    finding (sensory level) to its corresponding dermatome.
    """
    # In the clinical scenario, the key for localizing the spinal cord level
    # is the sensory deficit.
    sensory_level_landmark = "umbilicus"

    # A mapping of major anatomical landmarks to their spinal dermatome levels.
    dermatome_map = {
        "nipple line": "T4",
        "xiphoid process": "T6",
        "umbilicus": "T10",
        "inguinal region": "T12"
    }

    # The patient's presentation with ipsilateral motor/vibratory loss
    # and contralateral pain/temperature loss is classic Brown-Séquard syndrome.
    # The lesion level is determined by the upper limit of the sensory loss.
    lesion_level = dermatome_map.get(sensory_level_landmark)

    print("Step 1: The patient's symptoms point to Brown-Séquard syndrome (spinal cord hemisection).")
    print(f"Step 2: The sensory level for pain and temperature loss is at the '{sensory_level_landmark}'.")
    print(f"Step 3: The anatomical landmark '{sensory_level_landmark}' corresponds to the {lesion_level} dermatome.")
    print(f"Step 4: Therefore, the spinal cord injury is located at the {lesion_level} level.")

    # Match the determined level with the provided answer choices.
    answer_choices = {
        "A": "L4", "B": "L7", "C": "L5", "D": "T4", "E": "T6",
        "F": "None of the answer choices", "G": "T12", "H": "T10", "I": "C7"
    }

    final_answer_key = "Not Found"
    for key, value in answer_choices.items():
        if value == lesion_level:
            final_answer_key = key
            break

    print(f"\nConclusion: The injury at level {lesion_level} matches answer choice {final_answer_key}.")

find_injury_location()