def check_genetic_interactions():
    """
    Checks the genetic interactions based on the provided experimental data
    and evaluates the correctness of the final answer.
    """
    # --- Data from the question ---
    resistance = {
        "control": 100,
        "g1": 75,
        "g2": 0,
        "g3": 50,
        "g1g2": 0,
        "g1g3": 10,
        "g2g3": 0,
    }

    # The final answer provided is C. Let's find the text for option C from the candidate answers.
    # Most candidates agree on the text for the best-fit option.
    # e.g., Candidate 6, 12, 15: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    option_c_text = "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    
    # --- Analysis Functions ---

    def is_epistatic(gene_a, gene_b):
        """Checks if gene_a is epistatic to gene_b."""
        # Sort gene names to find the double mutant key, e.g., g1g2
        double_mutant_key = "".join(sorted((f"g{gene_a[-1]}", f"g{gene_b[-1]}")))
        
        # Epistasis condition: phenotype of double mutant == phenotype of one single mutant
        if resistance[double_mutant_key] == resistance[f"g{gene_a[-1]}"]:
            return True
        return False

    def check_tf_candidate(gene_name):
        """Checks if a gene is a likely transcription factor based on epistasis."""
        other_genes = [g for g in ['G1', 'G2', 'G3'] if g != gene_name]
        # A TF is likely epistatic to the genes it regulates.
        if is_epistatic(gene_name, other_genes[0]) and is_epistatic(gene_name, other_genes[1]):
            return True
        return False

    def check_g1_g3_interaction():
        """Analyzes the interaction between G1 and G3."""
        results = {}
        
        # Check for epistasis
        results['g1_epistatic_to_g3'] = is_epistatic('G1', 'G3')
        results['g3_epistatic_to_g1'] = is_epistatic('G3', 'G1')

        # Check for synergism (hallmark of gene redundancy)
        # Synergism: combined effect is greater than the sum of individual effects.
        loss_g1 = resistance['control'] - resistance['g1']  # 25
        loss_g3 = resistance['control'] - resistance['g3']  # 50
        loss_g1g3 = resistance['control'] - resistance['g1g3'] # 90
        
        # If the combined loss is greater than the sum of individual losses, it's synergistic.
        results['synergistic_interaction'] = loss_g1g3 > (loss_g1 + loss_g3)
        
        return results

    # --- Verification ---
    
    # Fact 1: Is G2 the transcription factor?
    is_g2_tf = check_tf_candidate('G2')
    if not is_g2_tf:
        return "Incorrect: The analysis that G2 is the transcription factor is wrong. G2 is not epistatic to both G1 and G3 according to the data."

    # Fact 2 & 3: Analyze G1-G3 interaction
    g1_g3_results = check_g1_g3_interaction()
    
    # Check the clauses of option C
    clause1_correct = is_g2_tf
    # "Gene redundancy" is inferred from the synergistic interaction
    clause2_correct = g1_g3_results['synergistic_interaction']
    clause3_correct = g1_g3_results['g1_epistatic_to_g3']

    # --- Final Evaluation ---
    
    # The provided answer's reasoning is that C is the "best fit" because the first two clauses are correct,
    # even though the third is incorrect. Let's verify this reasoning.
    
    if not clause1_correct:
        return "Incorrect: The first clause of option C, 'G2 is a transcription factor', is not supported by the data."
    
    if not clause2_correct:
        return "Incorrect: The second clause of option C, 'G1 and G3 show gene redundancy', is not supported. The interaction is not synergistic."
        
    if clause3_correct:
        return "Incorrect: The analysis states the third clause of option C is wrong, but the code finds it to be correct. There is a mismatch in reasoning."

    # If we reach here, it means the code confirms the provided answer's analysis:
    # Clause 1 is TRUE.
    # Clause 2 is TRUE.
    # Clause 3 is FALSE.
    # The decision to choose C as the best-fit option, despite its flawed third clause, is a valid judgment call because it captures the two most important conclusions from the experiment.
    
    return "Correct"

# Run the check
result = check_genetic_interactions()
print(result)