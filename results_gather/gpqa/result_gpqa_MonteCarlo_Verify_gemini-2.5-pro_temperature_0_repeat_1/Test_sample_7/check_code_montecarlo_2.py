def check_correctness():
    """
    This function checks the correctness of the LLM's answer by modeling the genetic interactions
    described in the question and verifying the claims made in the selected option and the reasoning.
    """
    # --- Data from the question ---
    # Resistance percentages for each genotype.
    resistance_data = {
        'wt': 100,
        'g1': 75,
        'g2': 0,
        'g3': 50,
        'g1g2': 0,
        'g1g3': 10,
        'g2g3': 0,
    }

    # --- Definitions of genetic concepts based on phenotypes ---

    def get_double_mutant_key(gene1, gene2):
        """Helper to create a sorted key for double mutants, e.g., ('g1', 'g2') -> 'g1g2'."""
        # Extracts numbers from gene names, sorts them, and reconstructs the key.
        genes = sorted([g.replace('g', '') for g in [gene1, gene2]])
        return f"g{genes[0]}g{genes[1]}"

    def is_epistatic(epistatic_gene, hypostatic_gene, data):
        """
        Checks if a gene is epistatic to another.
        This is true if the double mutant phenotype is the same as the epistatic gene's single mutant phenotype.
        """
        double_mutant_key = get_double_mutant_key(epistatic_gene, hypostatic_gene)
        return data[double_mutant_key] == data[epistatic_gene]

    def check_redundancy(gene1, gene2, data):
        """
        Checks for gene redundancy via a synergistic effect.
        This is true if the effect of the double mutant is greater than the sum of the effects of the single mutants.
        Effect is defined as the loss of resistance from wild-type (100%).
        """
        double_mutant_key = get_double_mutant_key(gene1, gene2)
        effect1 = 100 - data[gene1]
        effect2 = 100 - data[gene2]
        effect_double = 100 - data[double_mutant_key]
        return effect_double > (effect1 + effect2)

    def is_likely_tf(candidate_gene, other_genes, data):
        """
        Checks if a gene is a likely transcription factor.
        A key indicator in this context is being epistatic over other genes in the pathway,
        implying it acts upstream.
        """
        for other_gene in other_genes:
            if not is_epistatic(candidate_gene, other_gene, data):
                return False
        return True

    # --- Verification of the LLM's answer ---
    # The LLM chose option B: "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    # The LLM's reasoning is that the first two clauses are correct, and the third is incorrect, making B the best choice.
    # We will verify these three claims based on the data.

    # Claim 1: G2 is a transcription factor (because it's epistatic to G1 and G3).
    g2_is_tf = is_likely_tf('g2', ['g1', 'g3'], resistance_data)
    
    # Claim 2: G1 and G3 show gene redundancy.
    g1_g3_redundant = check_redundancy('g1', 'g3', resistance_data)

    # Claim 3: G1 is epistatic towards G3.
    g1_epistatic_to_g3 = is_epistatic('g1', 'g3', resistance_data)

    # The LLM's reasoning is sound if its analysis of the clauses in option B is correct.
    # The LLM states that clauses 1 and 2 are true, and clause 3 is false.
    if g2_is_tf and g1_g3_redundant and not g1_epistatic_to_g3:
        # The code's analysis matches the LLM's analysis.
        # The LLM correctly identified that G2 is a likely TF, G1 and G3 are redundant,
        # and that G1 is NOT epistatic to G3.
        # Based on this correct analysis, choosing B as the "best fit" is a valid conclusion.
        return "Correct"
    else:
        # If the code's analysis doesn't match the LLM's, the LLM's reasoning is flawed.
        reasons = []
        if g2_is_tf is not True:
            reasons.append("The answer's claim that G2 is a likely transcription factor (by being epistatic to others) is incorrect.")
        
        if g1_g3_redundant is not True:
            effect_g1 = 100 - resistance_data['g1']
            effect_g3 = 100 - resistance_data['g3']
            effect_g1g3 = 100 - resistance_data['g1g3']
            reasons.append(f"The answer's claim that G1 and G3 show redundancy is incorrect. The effect of g1g3 ({effect_g1g3}%) is not greater than the sum of individual effects ({effect_g1}% + {effect_g3}% = {effect_g1 + effect_g3}%).")

        if g1_epistatic_to_g3 is not False:
            reasons.append(f"The answer's claim that G1 is NOT epistatic to G3 is incorrect. The phenotype of g1g3 ({resistance_data['g1g3']}%) is equal to the phenotype of g1 ({resistance_data['g1']}%).")
        
        return "Incorrect. The reasoning provided in the answer is flawed for the following reason(s): " + " ".join(reasons)

# Execute the check and print the result.
print(check_correctness())