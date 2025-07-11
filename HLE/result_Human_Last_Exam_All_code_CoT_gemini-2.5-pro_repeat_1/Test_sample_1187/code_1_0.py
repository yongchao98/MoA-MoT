def solve_neurology_case():
    """
    Analyzes the clinical vignette to determine the location of the spinal cord injury.
    """
    # Patient's symptoms
    ipsilateral_motor_weakness = "Severe weakness in the right leg"
    ipsilateral_sensation_loss = "Complete loss of proprioceptive and vibratory sensation in the right leg"
    contralateral_sensation_loss = "Pain and temperature sensations are completely absent on her left side from the level of the umbilicus downward"

    # Explanation
    print("Step 1: Identify the syndrome based on the pattern of deficits.")
    print("The patient presents with a classic triad of symptoms:")
    print(f"- Ipsilateral motor loss (corticospinal tract): {ipsilateral_motor_weakness}")
    print(f"- Ipsilateral loss of proprioception/vibration (dorsal column): {ipsilateral_sensation_loss}")
    print(f"- Contralateral loss of pain/temperature (spinothalamic tract): {contralateral_sensation_loss}")
    print("This combination of findings is characteristic of Brown-Séquard syndrome, which is a hemisection of the spinal cord.\n")

    print("Step 2: Determine the anatomical level of the injury.")
    print("The sensory level for pain and temperature loss provides the most accurate localization of the lesion's vertical level.")
    print("The sensory loss is described as starting at the level of the umbilicus.")
    print("The dermatome corresponding to the umbilicus is T10.")
    print("Therefore, the spinal cord lesion is located at the T10 level.\n")

    print("Conclusion:")
    print("The patient's injury is at the T10 spinal cord level.")

solve_neurology_case()