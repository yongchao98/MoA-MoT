def check_answer_correctness():
    """
    Checks the correctness of the provided answer based on the experimental data.
    The chosen answer is B: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    """
    # Experimental data: resistance percentage
    resistance = {
        'WT': 100,
        'g1': 75,
        'g2': 0,
        'g3': 50,
        'g1g2': 0,
        'g1g3': 10,
        'g2g3': 0,
    }

    # --- Clause 1: G2 is a transcription factor ---
    # A TF is likely upstream, meaning it's epistatic to downstream genes.
    # Let's check if G2 is epistatic to G1 and G3.
    is_g2_epistatic_to_g1 = (resistance['g1g2'] == resistance['g2'])
    is_g2_epistatic_to_g3 = (resistance['g2g3'] == resistance['g2'])
    # The g2 phenotype is also the most severe (0%), consistent with a master regulator.
    is_g2_tf_correct = is_g2_epistatic_to_g1 and is_g2_epistatic_to_g3

    # --- Clause 2: G1 and G3 show gene redundancy ---
    # This is indicated by a synergistic interaction, where the combined loss is
    # greater than the sum of individual losses.
    loss_g1 = 100 - resistance['g1']  # 25
    loss_g3 = 100 - resistance['g3']  # 50
    loss_g1g3 = 100 - resistance['g1g3'] # 90
    is_redundancy_correct = loss_g1g3 > (loss_g1 + loss_g3)

    # --- Clause 3: G1 is epistatic towards G3 ---
    # This means the g1g3 phenotype should be the same as the g1 phenotype.
    is_g1_epistatic_g3_correct = (resistance['g1g3'] == resistance['g1'])

    # --- Final Verdict ---
    if is_g2_tf_correct and is_redundancy_correct and is_g1_epistatic_g3_correct:
        return "Correct"
    else:
        # The answer is incorrect because at least one clause is false.
        # Let's identify the specific reason.
        if not is_g1_epistatic_g3_correct:
            reason = (f"The statement 'G1 is epistatic towards G3' is false. "
                      f"For G1 to be epistatic, the phenotype of the g1g3 double mutant ({resistance['g1g3']}%) "
                      f"must be the same as the g1 single mutant ({resistance['g1']}%). Since they are different, "
                      f"the condition for epistasis is not met.")
            return f"Incorrect. The chosen answer B is not fully correct because {reason}"
        # Add other checks for completeness, though the above is the expected failure.
        elif not is_g2_tf_correct:
            return "Incorrect. The claim 'G2 is a transcription factor' is not supported by the data."
        elif not is_redundancy_correct:
            return "Incorrect. The claim 'G1 and G3 show gene redundancy' is not supported by the data."

# Run the check
result = check_answer_correctness()
print(result)