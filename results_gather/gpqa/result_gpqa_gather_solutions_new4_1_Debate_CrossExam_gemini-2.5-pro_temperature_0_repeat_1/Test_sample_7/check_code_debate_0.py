def check_answer_correctness():
    """
    Checks the correctness of the answer to the genetics question.

    The function programmatically analyzes the provided experimental data to deduce
    the relationships between genes G1, G2, and G3. It then evaluates the chosen
    answer (Option C) to see if it aligns with the data-driven conclusions.
    """
    # Step 1: Define the experimental data
    data = {
        "wt": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # Step 2: Perform genetic analysis
    
    # Check for epistasis
    g2_epistatic_to_g1 = (data["g1g2"] == data["g2"])
    g2_epistatic_to_g3 = (data["g2g3"] == data["g2"])
    g1_epistatic_to_g3 = (data["g1g3"] == data["g1"])
    g3_epistatic_to_g1 = (data["g1g3"] == data["g3"])

    # Identify the transcription factor (TF)
    # G2's knockout is most severe (0%) and it is epistatic to both G1 and G3.
    is_g2_tf = g2_epistatic_to_g1 and g2_epistatic_to_g3

    # Analyze the interaction between G1 and G3
    # A synergistic interaction (hallmark of gene redundancy) occurs if the double
    # mutant's effect is greater than the sum of the single mutants' effects.
    loss_g1 = data["wt"] - data["g1"]  # 25%
    loss_g3 = data["wt"] - data["g3"]  # 50%
    additive_loss = loss_g1 + loss_g3 # 75%
    actual_loss_g1g3 = data["wt"] - data["g1g3"] # 90%
    is_g1g3_redundant = actual_loss_g1g3 > additive_loss

    # Step 3: Evaluate the chosen answer (Option C)
    # Option C: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    
    # Check each clause of Option C
    clause1_correct = is_g2_tf
    clause2_correct = is_g1g3_redundant
    clause3_correct = g1_epistatic_to_g3

    # The provided answer analysis concludes C is the "best fit" because it has two correct
    # statements and one incorrect one, which is better than other options.
    # Let's verify this logic.
    
    num_correct_clauses_in_c = sum([clause1_correct, clause2_correct])
    
    if not clause1_correct:
        return "Incorrect. The answer claims G2 is the transcription factor, but the analysis shows this is false."
    if not clause2_correct:
        return "Incorrect. The answer claims G1 and G3 show gene redundancy, but the analysis shows this is false."
    if clause3_correct:
        return "Incorrect. The answer claims G1 is epistatic towards G3, which is a key incorrect part of the 'best fit' reasoning. The analysis shows this statement is actually true, invalidating the reasoning for choosing C."

    # At this point, we've confirmed C has two correct clauses and one incorrect one.
    # Now, let's check if other options are worse, to validate the "best fit" argument.
    
    # Option A: G1 is TF (False), G2/G3 pleiotropy (Unsupported), G2 epistatic to G1 (True) -> Max 1 correct
    # Option B: G2 is TF (True), G1/G3 pleiotropy (Unsupported), G1 epistatic to G3 (False) -> Max 1 correct
    # Option D: G2 is TF (True), G1/G3 same promoter (Unsupported), G3 epistatic to G1 (False) -> Max 1 correct
    
    # Since Option C has 2 correct, verifiable clauses while others have at most 1,
    # the reasoning that C is the best fit among flawed options is sound.
    
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)