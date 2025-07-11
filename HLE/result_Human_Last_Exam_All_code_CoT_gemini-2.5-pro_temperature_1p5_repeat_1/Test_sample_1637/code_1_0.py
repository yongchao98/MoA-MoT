def calculate_solubility_score():
    """
    Calculates and compares hydrophobicity scores for the original and
    a modified chemical probe to predict solubility changes.

    This is a simplified model where:
    - Positive scores represent hydrophobic (water-repelling) character.
    - Negative scores represent hydrophilic (water-attracting) character.
    - A lower total score suggests better water solubility.
    """

    # Assigning arbitrary hydrophobicity scores to molecular fragments
    thioxanthenone_core_score = 7.0
    chlorohexyl_chain_score = 4.0
    # The original linker has some hydrophilic ether and amide groups
    original_linker_score = -5.0

    # The proposed PEG linker is significantly more hydrophilic
    # This represents replacing the original linker with a PEG-based one
    peg_linker_score = -9.0

    # Calculate the total score for the original probe
    original_probe_total_score = (thioxanthenone_core_score +
                                  chlorohexyl_chain_score +
                                  original_linker_score)

    # Calculate the total score for the new probe with a PEG linker
    new_probe_total_score = (thioxanthenone_core_score +
                             chlorohexyl_chain_score +
                             peg_linker_score)

    print("--- Predicting Probe Solubility ---")
    print("A lower score indicates higher predicted water solubility.\n")

    # --- Output for the Original Probe ---
    print("Original Probe Hydrophobicity Equation:")
    # Printing each number in the equation
    print(f"Score = {thioxanthenone_core_score} (Thioxanthenone) + "
          f"{chlorohexyl_chain_score} (Chlorohexyl) + "
          f"{original_linker_score} (Amide Linker)")
    print(f"Original Probe Total Score = {original_probe_total_score:.1f}\n")


    # --- Output for the New Probe ---
    print("Proposed Probe with PEG Linker Hydrophobicity Equation:")
    # Printing each number in the equation
    print(f"Score = {thioxanthenone_core_score} (Thioxanthenone) + "
          f"{chlorohexyl_chain_score} (Chlorohexyl) + "
          f"{peg_linker_score} (PEG Linker)")
    print(f"New Probe Total Score = {new_probe_total_score:.1f}\n")

    # --- Conclusion ---
    print("Conclusion:")
    if new_probe_total_score < original_probe_total_score:
        print("The new probe with the PEG linker has a lower hydrophobicity score.")
        print("This suggests it will have significantly better water solubility, "
              "which should solve the precipitation problem.")
    else:
        print("Based on this model, the change might not improve solubility.")

# Run the calculation and print the results
calculate_solubility_score()