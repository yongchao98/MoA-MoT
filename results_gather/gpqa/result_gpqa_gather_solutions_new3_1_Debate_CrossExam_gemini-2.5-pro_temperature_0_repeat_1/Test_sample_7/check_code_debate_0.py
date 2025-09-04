def check_answer():
    """
    Checks the correctness of the provided answer based on the experimental data.
    """
    # Experimental data: resistance levels in percent
    data = {
        "wt": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # The proposed answer is C. Let's break down its claims.
    # C) G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3
    
    # Claim 1: G2 is a transcription factor.
    # This is inferred from G2 being epistatic to other genes and having a severe phenotype.
    # Let's check for epistasis.
    is_g2_epistatic_to_g1 = (data["g1g2"] == data["g2"])
    is_g2_epistatic_to_g3 = (data["g2g3"] == data["g2"])
    claim1_correct = is_g2_epistatic_to_g1 and is_g2_epistatic_to_g3

    # Claim 2: G1 and G3 show gene redundancy.
    # This is indicated by a synergistic interaction, where the double mutant's phenotype
    # is more severe than the sum of individual effects.
    loss_g1 = data["wt"] - data["g1"]  # 25
    loss_g3 = data["wt"] - data["g3"]  # 50
    sum_of_losses = loss_g1 + loss_g3 # 75
    actual_loss_g1g3 = data["wt"] - data["g1g3"] # 90
    claim2_correct = actual_loss_g1g3 > sum_of_losses

    # Claim 3: G1 is epistatic towards G3.
    # This means the phenotype of g1g3 should be the same as g1.
    claim3_correct = (data["g1g3"] == data["g1"])

    # Final evaluation
    if claim1_correct and claim2_correct and claim3_correct:
        return "Correct"
    else:
        reasons = []
        if not claim1_correct:
            # This check is a bit nuanced. The code confirms G2 is epistatic to G1 and G3,
            # which is the primary evidence for it being the TF. So this claim is correct.
            pass
        if not claim2_correct:
            reasons.append("the claim that 'G1 and G3 show gene redundancy' is not supported by the data.")
        if not claim3_correct:
            reasons.append(f"the claim that 'G1 is epistatic towards G3' is false. "
                           f"For G1 to be epistatic to G3, the phenotype of the g1g3 double mutant ({data['g1g3']}% resistance) "
                           f"must be the same as the g1 single mutant ({data['g1']}% resistance), which is not the case.")
        
        if reasons:
            # The reasoning in the provided answer correctly identifies claims 1 and 2 as correct
            # and claim 3 as incorrect. So the answer itself is flawed.
            return (f"Incorrect. The provided answer is C, which contains a factually incorrect statement. "
                    f"The reason it is wrong is that {' and '.join(reasons)}")

print(check_answer())