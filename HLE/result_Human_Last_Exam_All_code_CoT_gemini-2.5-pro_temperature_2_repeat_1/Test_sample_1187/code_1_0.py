def diagnose_spinal_injury():
    """
    Analyzes clinical signs to determine the location of a spinal cord injury.
    """
    
    # Patient's key symptoms
    right_leg_weakness = "Ipsilateral motor deficit (Corticospinal tract)"
    right_leg_proprioception_vibration_loss = "Ipsilateral proprioception/vibration loss (Dorsal columns)"
    left_side_pain_temp_loss = "Contralateral pain/temperature loss (Spinothalamic tract)"
    sensory_level = "Umbilicus"

    # Step 1: Identify the syndrome
    print("Clinical Analysis:")
    print("1. The patient has weakness in the right leg (ipsilateral to the right-sided stab wound). This points to damage to the right corticospinal tract.")
    print("2. The patient has loss of proprioception and vibration in the right leg (ipsilateral). This indicates damage to the right dorsal columns.")
    print("3. The patient has loss of pain and temperature sensation on the left side (contralateral). This indicates damage to the right spinothalamic tract, as its fibers cross over shortly after entering the cord.")
    print("\nConclusion 1: This classic triad of symptoms is known as Brown-Séquard Syndrome (spinal cord hemisection), in this case, on the right side.")

    # Step 2: Determine the neurological level
    # Dermatome map for key landmarks
    dermatomes = {
        "Nipple line": "T4",
        "Xiphoid process": "T6",
        "Umbilicus": "T10",
        "Inguinal region": "T12"
    }
    
    lesion_level_dermatome = dermatomes[sensory_level]
    
    print("\nLocalization of the Lesion:")
    print(f"The sensory level for pain and temperature loss is noted at the {sensory_level}.")
    print(f"The dermatome corresponding to the umbilicus is T10.")
    print("Since the spinothalamic tract is affected, this sensory level directly corresponds to the level of the spinal cord injury.")
    
    # Step 3: Match with answer choices
    print(f"\nConclusion 2: The injury is located at the T10 spinal level.")

diagnose_spinal_injury()