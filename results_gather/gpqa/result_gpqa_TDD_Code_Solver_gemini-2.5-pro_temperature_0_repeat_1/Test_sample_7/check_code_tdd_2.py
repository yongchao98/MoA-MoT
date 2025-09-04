def check_correctness_of_answer():
    """
    This function checks the correctness of the provided LLM answer by programmatically
    verifying the genetic interactions based on the experimental data.
    """
    # --- Data from the question ---
    resistance = {
        'wt': 100,
        'g1': 75,
        'g2': 0,
        'g3': 50,
        'g1g2': 0,
        'g1g3': 10,
        'g2g3': 0
    }

    # --- List to store reasons for incorrectness ---
    error_messages = []

    # --- Part 1: Verify the Transcription Factor (TF) identification ---
    # The TF should be epistatic to the other genes.
    # A gene (e.g., gA) is epistatic to another (gB) if the phenotype of the double mutant (gAgB)
    # is the same as the phenotype of the single epistatic mutant (gA).
    is_g2_epistatic_to_g1 = (resistance['g1g2'] == resistance['g2'])
    is_g2_epistatic_to_g3 = (resistance['g2g3'] == resistance['g2'])
    
    # The answer claims G2 is the TF. Let's check if this holds.
    if not (is_g2_epistatic_to_g1 and is_g2_epistatic_to_g3):
        error_messages.append(
            "Constraint Violated: The answer identifies G2 as the transcription factor due to epistasis. "
            f"However, the data does not support G2 being epistatic to both G1 and G3. "
            f"G2 epistatic to G1: {is_g2_epistatic_to_g1}. G2 epistatic to G3: {is_g2_epistatic_to_g3}."
        )
    
    # Let's also check if other genes could be the TF, which would invalidate the answer's conclusion.
    is_g1_epistatic_to_g2 = (resistance['g1g2'] == resistance['g1']) # 0 == 75 -> False
    if is_g1_epistatic_to_g2:
        error_messages.append(
            "Constraint Violated: The answer claims G2 is the TF, but the data shows G1 is epistatic to G2, "
            "which contradicts the proposed genetic pathway."
        )

    # --- Part 2: Verify the interaction between G1 and G3 ---
    # The answer claims gene redundancy, which is indicated by a synergistic interaction.
    # Synergy: The combined effect is greater than the sum of individual effects.
    loss_g1 = resistance['wt'] - resistance['g1']
    loss_g3 = resistance['wt'] - resistance['g3']
    expected_additive_loss = loss_g1 + loss_g3
    observed_double_mutant_loss = resistance['wt'] - resistance['g1g3']

    if not (observed_double_mutant_loss > expected_additive_loss):
        error_messages.append(
            "Constraint Violated: The answer claims G1 and G3 show gene redundancy based on a synergistic effect. "
            f"However, the observed loss in the double mutant ({observed_double_mutant_loss}%) is not greater than "
            f"the expected additive loss ({expected_additive_loss}%)."
        )

    # --- Part 3: Evaluate the specific claims in the chosen answer (C) ---
    # C) G2 is a TF, G1/G3 show gene redundancy, G1 is epistatic towards G3.
    # The LLM's answer correctly notes that the "G1 is epistatic" part is false. Let's confirm.
    is_g1_epistatic_to_g3 = (resistance['g1g3'] == resistance['g1']) # 10 == 75 -> False
    if is_g1_epistatic_to_g3:
        error_messages.append(
            "Reasoning Flaw: The provided answer correctly states that the 'G1 is epistatic towards G3' clause in option C is false. "
            "However, this check shows that the clause is actually true, meaning the answer's reasoning for selecting C is flawed."
        )

    # --- Final Verdict ---
    if not error_messages:
        # The code confirms that the logical steps and conclusions in the provided answer are correct.
        # 1. G2 is correctly identified as the TF due to epistasis.
        # 2. The G1-G3 interaction is correctly identified as synergistic (gene redundancy).
        # 3. The answer correctly identifies that the third clause of option C is false, but that C is still the best option.
        return "Correct"
    else:
        return "Incorrect. " + " ".join(error_messages)

# Run the check
result = check_correctness_of_answer()
print(result)