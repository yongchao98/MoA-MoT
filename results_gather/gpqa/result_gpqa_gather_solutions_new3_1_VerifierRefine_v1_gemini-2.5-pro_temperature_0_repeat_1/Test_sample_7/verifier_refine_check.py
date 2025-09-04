import re

def check_answer_correctness():
    """
    This function checks the correctness of the final answer (Option D) based on the provided genetic data.
    The function evaluates each clause of the answer statement.
    """
    # 1. Store the experimental data from the question
    data = {
        'WT': 100,
        'g1': 75,
        'g2': 0,
        'g3': 50,
        'g1g2': 0,
        'g1g3': 10,
        'g2g3': 0
    }

    # 2. Define helper functions to test genetic principles based on the data

    def is_epistatic(epistatic_gene, masked_gene):
        """
        Checks if epistatic_gene is epistatic to masked_gene.
        The phenotype of the double mutant must match the phenotype of the epistatic gene's single mutant.
        e.g., is_epistatic('G2', 'G1') checks if phenotype(g1g2) == phenotype(g2)
        """
        epistatic_mutant_key = 'g' + epistatic_gene[1:]
        masked_mutant_key = 'g' + masked_gene[1:]
        
        # Create possible keys for the double mutant, e.g., 'g1g2' or 'g2g1'
        double_mutant_key1 = epistatic_mutant_key + masked_mutant_key[1:]
        double_mutant_key2 = masked_mutant_key + epistatic_mutant_key[1:]
        
        double_mutant_phenotype = data.get(double_mutant_key1, data.get(double_mutant_key2))
        epistatic_phenotype = data.get(epistatic_mutant_key)

        if double_mutant_phenotype is None or epistatic_phenotype is None:
            return False

        return double_mutant_phenotype == epistatic_phenotype

    def check_transcription_factor(gene, other_genes):
        """
        Checks if a gene is a likely transcription factor (TF).
        Based on the problem, a TF is upstream (epistatic) and its knockout has a severe phenotype.
        """
        # Check for severe phenotype (complete loss of function)
        mutant_key = 'g' + gene[1:]
        if data.get(mutant_key) != 0:
            return False
            
        # Check if it's epistatic to all other specified genes
        is_upstream = all(is_epistatic(gene, other) for other in other_genes)
        
        return is_upstream

    def check_gene_redundancy(gene_A, gene_B):
        """
        Checks for gene redundancy via synergistic interaction.
        Synergy is defined as the loss of function in the double mutant being greater
        than the sum of the losses in the single mutants.
        """
        loss_A = 100 - data['g' + gene_A[1:]]
        loss_B = 100 - data['g' + gene_B[1:]]
        
        double_mutant_key1 = 'g' + gene_A[1:] + 'g' + gene_B[1:]
        double_mutant_key2 = 'g' + gene_B[1:] + 'g' + gene_A[1:]
        double_mutant_phenotype = data.get(double_mutant_key1, data.get(double_mutant_key2))
        loss_double = 100 - double_mutant_phenotype

        return loss_double > (loss_A + loss_B)

    # 3. Evaluate the claims in the provided answer (Option D)
    # The answer to check is: "D) G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3"
    
    errors = []

    # Claim 1: "G2 is a transcription factor"
    claim1_correct = check_transcription_factor('G2', ['G1', 'G3'])
    if not claim1_correct:
        errors.append("The claim 'G2 is a transcription factor' is not supported.")

    # Claim 2: "G1 and G3 show gene redundancy"
    claim2_correct = check_gene_redundancy('G1', 'G3')
    if not claim2_correct:
        errors.append("The claim 'G1 and G3 show gene redundancy' is false.")

    # Claim 3: "G1 is epistatic towards G3"
    claim3_correct = is_epistatic('G1', 'G3')
    if not claim3_correct:
        errors.append(f"The claim 'G1 is epistatic towards G3' is false. For G1 to be epistatic to G3, the phenotype of the double mutant g1g3 ({data['g1g3']}%) must be the same as the phenotype of the single mutant g1 ({data['g1']}%). Since {data['g1g3']}% is not equal to {data['g1']}%, this condition is not met.")

    # 4. Formulate the final output
    if not errors:
        return "Correct"
    else:
        # The provided answer is a compound statement. If any part is false, the whole statement is false.
        return "Incorrect. The provided answer (Option D) contains a factually incorrect statement: " + " ".join(errors)

# Execute the check and print the result.
result = check_answer_correctness()
print(result)