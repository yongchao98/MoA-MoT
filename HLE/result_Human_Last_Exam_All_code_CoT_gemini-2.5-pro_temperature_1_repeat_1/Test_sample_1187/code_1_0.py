def explain_diagnosis():
    """
    Explains the reasoning for identifying the location of the spinal cord injury.
    """
    
    # Patient's symptoms
    right_leg_weakness = "Ipsilateral motor loss (Corticospinal Tract)"
    right_leg_sensation_loss = "Ipsilateral loss of proprioception/vibration (Dorsal Columns)"
    left_side_pain_loss = "Contralateral loss of pain/temperature (Spinothalamic Tract)"
    sensory_level = "Umbilicus"
    
    # Diagnosis
    syndrome = "Brown-Séquard Syndrome (Spinal Cord Hemisection)"
    
    # Dermatome Level
    dermatome_level = "T10"
    
    print("Clinical Analysis Steps:")
    print("1. The patient's presentation with ipsilateral motor weakness, ipsilateral loss of proprioception/vibration, and contralateral loss of pain/temperature is characteristic of Brown-Séquard Syndrome.")
    print(f"2. The key to locating the lesion vertically is the sensory level. The loss of pain and temperature begins at the umbilicus.")
    print(f"3. The dermatome that supplies the skin at the level of the umbilicus is T10.")
    print(f"4. Therefore, the injury is located at the T10 spinal cord level.")
    print("\nConclusion:")
    print(f"The location of the patient's injury is at the {dermatome_level} spinal segment.")

explain_diagnosis()