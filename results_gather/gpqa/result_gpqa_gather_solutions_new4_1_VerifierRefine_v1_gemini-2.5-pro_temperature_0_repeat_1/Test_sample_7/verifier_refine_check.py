import re

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer for the genetics question.
    It does this by:
    1. Storing the experimental data.
    2. Defining functions to test genetic principles (epistasis, synergy/redundancy, TF identification).
    3. Parsing each multiple-choice option into its constituent clauses.
    4. Evaluating the correctness of each clause based on the data.
    5. Scoring each option to find the "best" fit, even if imperfect.
    6. Comparing the determined best option with the LLM's answer.
    """
    # --- Data from the question and the LLM's proposed answer ---
    resistance = {
        'WT': 100,
        'g1': 75,
        'g2': 0,
        'g3': 50,
        'g1g2': 0,
        'g1g3': 10,
        'g2g3': 0
    }
    llm_answer = 'B'

    # --- Helper functions to check genetic principles ---
    def is_epistatic(epistatic_gene, masked_gene, data):
        """Checks if epistatic_gene masks masked_gene.
        The double mutant phenotype must match the epistatic gene's single mutant phenotype.
        """
        g_epistatic_name = f"g{epistatic_gene[-1]}"
        g_masked_name = f"g{masked_gene[-1]}"
        
        # Ensure correct key for double mutant (e.g., g1g2, not g2g1)
        key_parts = sorted([int(g_epistatic_name[1:]), int(g_masked_name[1:])])
        double_mutant_key = f"g{key_parts[0]}g{key_parts[1]}"
        
        return data[double_mutant_key] == data[g_epistatic_name]

    def shows_synergy_for_redundancy(gene1, gene2, data):
        """Checks for synergistic interaction, a hallmark of gene redundancy.
        The combined loss of function must be greater than the sum of individual losses.
        """
        g1_name = f"g{gene1[-1]}"
        g2_name = f"g{gene2[-1]}"
        
        loss_g1 = data['WT'] - data[g1_name]
        loss_g2 = data['WT'] - data[g2_name]
        sum_of_individual_losses = loss_g1 + loss_g2
        
        key_parts = sorted([int(g1_name[1:]), int(g2_name[1:])])
        double_mutant_key = f"g{key_parts[0]}g{key_parts[1]}"
        observed_loss = data['WT'] - data[double_mutant_key]
        
        return observed_loss > sum_of_individual_losses

    def is_likely_tf(gene, other_genes, data):
        """A likely TF has the most severe phenotype (0%) and is epistatic to others."""
        gene_name = f"g{gene[-1]}"
        # Condition 1: Most severe phenotype (complete loss of function)
        if data[gene_name] != 0:
            return False
        # Condition 2: Epistatic to other genes in the pathway
        for other in other_genes:
            if not is_epistatic(gene, other, data):
                return False
        return True

    # --- Analysis of the options ---
    options = {
        'A': "G1 is a transcription factor, G2 and G3 show pleiotropy, G2 is epistatic towards G1",
        'B': "G2 is a transcription factor, G1 and G3 show gene redundancy, G1 is epistatic towards G3",
        'C': "G2 is a transcription factor, G1 and G3 show pleiotropy, G1 is epistatic towards G3",
        'D': "G2 is a transcription factor, G1 and G3 has the same promoter, G3 is epistatic towards G1"
    }
    
    analysis_results = {}
    
    for option_key, option_text in options.items():
        clauses = [c.strip() for c in option_text.split(',')]
        results = []
        
        # Clause 1: Transcription Factor
        tf_match = re.search(r"(G\d) is a transcription factor", clauses[0])
        gene = tf_match.group(1)
        other_genes = [g for g in ['G1', 'G2', 'G3'] if g != gene]
        results.append({'is_correct': is_likely_tf(gene, other_genes, resistance)})

        # Clause 2: G1/G3 interaction
        interaction_clause = clauses[1]
        if "gene redundancy" in interaction_clause:
            results.append({'is_correct': shows_synergy_for_redundancy('G1', 'G3', resistance)})
        elif "pleiotropy" in interaction_clause or "same promoter" in interaction_clause:
            # These cannot be determined from the given data
            results.append({'is_correct': 'Unverifiable'})
        else:
            results.append({'is_correct': False})

        # Clause 3: Epistasis
        epistasis_match = re.search(r"(G\d) is epistatic towards (G\d)", clauses[2])
        epistatic_gene = epistasis_match.group(1)
        masked_gene = epistasis_match.group(2)
        results.append({'is_correct': is_epistatic(epistatic_gene, masked_gene, resistance)})
            
        analysis_results[option_key] = results

    # --- Determine the best option based on analysis ---
    best_option = None
    max_score = -1
    
    for option_key, results in analysis_results.items():
        score = 0
        # Disqualify if the main premise (TF identification) is wrong. This is a critical error.
        if results[0]['is_correct'] == False:
            continue
        
        # Score based on correct, verifiable statements
        for res in results:
            if res['is_correct'] == True:
                score += 1
        
        if score > max_score:
            max_score = score
            best_option = option_key

    # --- Final check ---
    if llm_answer == best_option:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect. My analysis shows that '{best_option}' is the best choice.\n\n"
        reason += "Here is a breakdown of the correct conclusions from the data:\n"
        reason += "1. G2 is the transcription factor: The g2 mutant has 0% resistance and is epistatic to both G1 and G3.\n"
        reason += "2. G1 and G3 show gene redundancy: The g1g3 double mutant (10% resistance) has a much more severe phenotype than expected from adding the single mutant effects. This synergistic interaction is a hallmark of redundancy.\n"
        reason += "3. There is NO epistasis between G1 and G3: The g1g3 phenotype (10%) is not the same as g1 (75%) or g3 (50%).\n\n"
        reason += "Evaluating the options based on these facts:\n"
        reason += "- Option A: Incorrectly identifies G1 as the TF.\n"
        reason += "- Option B: Correctly identifies G2 as TF (Correct) and G1/G3 as redundant (Correct). The final clause is wrong. (Score: 2)\n"
        reason += "- Option C: Correctly identifies G2 as TF (Correct), but 'pleiotropy' is a less accurate/verifiable description than 'redundancy', and the final clause is wrong. (Score: 1)\n"
        reason += "- Option D: Correctly identifies G2 as TF (Correct), but 'same promoter' is unverifiable, and the final clause is wrong. (Score: 1)\n\n"
        reason += f"Conclusion: Option B is the best fit as it contains the two most significant and correct conclusions. The provided answer '{llm_answer}' does not align with this analysis."
        return reason

# The code will return "Correct" because its analysis aligns with the provided answer's logic:
# Option B correctly identifies the two most important relationships (G2 as the upstream TF and the redundant/synergistic nature of G1 and G3), making it the best choice despite the final clause being incorrect.
print(check_correctness())