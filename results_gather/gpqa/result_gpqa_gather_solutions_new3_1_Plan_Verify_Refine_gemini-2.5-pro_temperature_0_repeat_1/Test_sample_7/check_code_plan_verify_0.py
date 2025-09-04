def check_genetic_conclusions():
    """
    Checks the correctness of the proposed answer based on the experimental data.
    The proposed answer is A: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    """
    # Store the experimental data in a dictionary
    resistance = {
        "control": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # --- Check Claim 1: G2 is a transcription factor (inferred from being an upstream epistatic regulator) ---
    # G2 is epistatic to G1 if the g1g2 phenotype matches the g2 phenotype.
    is_g2_epistatic_to_g1 = (resistance["g1g2"] == resistance["g2"])
    # G2 is epistatic to G3 if the g2g3 phenotype matches the g2 phenotype.
    is_g2_epistatic_to_g3 = (resistance["g2g3"] == resistance["g2"])
    # The knockout of an upstream TF usually has the most severe phenotype.
    g2_has_most_severe_phenotype = resistance["g2"] <= resistance["g1"] and resistance["g2"] <= resistance["g3"]
    
    claim1_supported = is_g2_epistatic_to_g1 and is_g2_epistatic_to_g3 and g2_has_most_severe_phenotype

    # --- Check Claim 2: G1 and G3 show gene redundancy (inferred from synergistic interaction) ---
    # Synergy: The combined loss of function is greater than the sum of individual losses.
    loss_from_g1 = 100 - resistance["g1"]  # 25%
    loss_from_g3 = 100 - resistance["g3"]  # 50%
    additive_loss = loss_from_g1 + loss_from_g3 # 75%
    observed_loss_in_g1g3 = 100 - resistance["g1g3"] # 90%
    
    claim2_supported = observed_loss_in_g1g3 > additive_loss

    # --- Check Claim 3: G1 is epistatic towards G3 ---
    # This requires the g1g3 phenotype to be identical to the g1 phenotype.
    claim3_supported = (resistance["g1g3"] == resistance["g1"])

    # --- Final Verdict ---
    # The answer is only correct if all its claims are true.
    if not claim1_supported:
        return "Incorrect. The claim 'G2 is a transcription factor' is not supported by the epistasis data."
    
    if not claim2_supported:
        return f"Incorrect. The claim 'G1 and G3 show gene redundancy' is not supported. The observed loss in the g1g3 double mutant ({observed_loss_in_g1g3}%) is not greater than the additive loss from single mutants ({additive_loss}%)."

    if not claim3_supported:
        return f"Incorrect. The claim 'G1 is epistatic towards G3' is not satisfied. For G1 to be epistatic to G3, the resistance of the g1g3 double mutant must be the same as the g1 single mutant. However, the resistance of g1g3 is {resistance['g1g3']}% while the resistance of g1 is {resistance['g1']}%."

    return "Correct"

# Run the check and print the result
result = check_genetic_conclusions()
print(result)