def check_genetics_answer():
    """
    Checks the correctness of the proposed answer based on the experimental data.
    """
    # Experimental data: resistance percentage
    data = {
        "WT": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # The final answer from the LLM to be checked is C.
    # Option C claims: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    
    # --- Verification of Claim 1: G2 is a transcription factor ---
    # A TF is likely upstream and epistatic to downstream genes.
    # Check if G2 is epistatic to G1
    g2_epistatic_to_g1 = (data["g1g2"] == data["g2"])
    # Check if G2 is epistatic to G3
    g2_epistatic_to_g3 = (data["g2g3"] == data["g2"])
    # Check if g2 mutant has the most severe phenotype
    g2_most_severe = data["g2"] <= data["g1"] and data["g2"] <= data["g3"]
    
    claim1_correct = g2_epistatic_to_g1 and g2_epistatic_to_g3 and g2_most_severe

    # --- Verification of Claim 2: G1 and G3 show gene redundancy ---
    # Gene redundancy is indicated by a synergistic interaction, where the double mutant's
    # loss of function is greater than the sum of the single mutants' losses of function.
    loss_g1 = 100 - data["g1"]  # 25%
    loss_g3 = 100 - data["g3"]  # 50%
    sum_of_losses = loss_g1 + loss_g3 # 75%
    
    loss_g1g3 = 100 - data["g1g3"] # 90%
    
    claim2_correct = loss_g1g3 > sum_of_losses

    # --- Verification of Claim 3: G1 is epistatic towards G3 ---
    # For G1 to be epistatic to G3, the g1g3 phenotype must equal the g1 phenotype.
    claim3_correct = (data["g1g3"] == data["g1"])

    # --- Final Verdict ---
    if claim1_correct and claim2_correct and claim3_correct:
        return "Correct"
    else:
        reasons = []
        if not claim1_correct:
            # This claim is actually correct, but we include the check for completeness.
            reasons.append("Claim 'G2 is a transcription factor' is not fully supported.")
        if not claim2_correct:
            reasons.append("Claim 'G1 and G3 show gene redundancy' is not supported. The interaction is not synergistic.")
        if not claim3_correct:
            reasons.append(f"Claim 'G1 is epistatic towards G3' is not satisfied. The phenotype of the g1g3 double mutant ({data['g1g3']}% resistance) is not the same as the phenotype of the g1 single mutant ({data['g1']}% resistance).")
        
        # Since the overall answer C contains a false statement, it is incorrect.
        # We return the most direct and verifiable reason for its incorrectness.
        return f"The answer is incorrect because not all of its statements are true. Specifically, the following condition is not satisfied:\n- {reasons[-1]}"

# Run the check
result = check_genetics_answer()
print(result)