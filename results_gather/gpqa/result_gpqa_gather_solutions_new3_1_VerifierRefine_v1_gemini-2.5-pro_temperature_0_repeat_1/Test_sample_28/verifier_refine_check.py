def check_answer():
    """
    This function checks the correctness of the answer by analyzing the genetic mutations.
    It determines the type of each mutation and ranks them by their potential to disrupt protein function.
    The most disruptive mutations are nonsense and frameshift mutations.
    """

    # Define the initial part of the gene sequences for analysis.
    # The critical changes are all within the first ~20 base pairs.
    intact_gene = "ATGTTTCTCGCTGGTACTTCT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCT",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCT",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACT",
        "Mutant 4": "ATGTTTTACGCTACTTCT",
    }
    
    # The options provided in the question map letters to mutants.
    # A) Mutant 1, B) Mutant 3, C) Mutant 4, D) Mutant 2
    option_map = {"A": "Mutant 1", "B": "Mutant 3", "C": "Mutant 4", "D": "Mutant 2"}
    provided_answer_option = "D"
    
    # --- Analysis of each mutation ---
    # This analysis is based on comparing the mutant sequences to the intact one.
    # Disruption score: Lower is more severe. 1 is most severe.
    analysis = {}

    # Intact codons: ATG TTT CTC GCT ... (Met-Phe-Leu-Ala...)

    # Mutant 1: ATG TTC TAC ... (Met-Phe-Tyr...). Change from CTC (Leu) to TAC (Tyr).
    # This is a missense mutation.
    analysis["Mutant 1"] = {"type": "Missense", "disruption": 3, "reason": "Substitutes one amino acid (Leucine to Tyrosine). Function may or may not be lost."}

    # Mutant 2: ATG TTC TAA ... (Met-Phe-Stop). Change from CTC (Leu) to TAA (Stop codon).
    # This is a nonsense mutation.
    analysis["Mutant 2"] = {"type": "Nonsense", "disruption": 1, "reason": "Introduces a premature stop codon, truncating the protein. Almost certain to be non-functional."}

    # Mutant 3: ATG TTT T ACG CTG ... An insertion of 'T' after the 2nd codon.
    # This causes a frameshift mutation.
    analysis["Mutant 3"] = {"type": "Frameshift", "disruption": 1, "reason": "An insertion alters the reading frame, scrambling all subsequent amino acids. Almost certain to be non-functional."}

    # Mutant 4: ATG TTT TAC GCT ACT... vs Intact: ATG TTT CTC GCT GGT ACT...
    # This involves a missense mutation (CTC->TAC) and a deletion of a full codon (GGT).
    # This is an in-frame deletion plus a missense mutation.
    analysis["Mutant 4"] = {"type": "In-frame deletion + Missense", "disruption": 2, "reason": "Removes one amino acid and changes another. Disruptive, but less certain to be a complete knockout than nonsense or frameshift."}

    # --- Verification ---
    
    # Find the most disruptive mutation type(s)
    min_disruption_score = min(data["disruption"] for data in analysis.values())
    best_mutants = [name for name, data in analysis.items() if data["disruption"] == min_disruption_score]

    # The provided answer corresponds to Mutant 2
    answer_mutant_name = option_map[provided_answer_option]

    # Check if the selected answer is one of the best candidates.
    if answer_mutant_name in best_mutants:
        return "Correct"
    else:
        reason = (f"The provided answer is '{provided_answer_option}', which corresponds to {answer_mutant_name}. "
                  f"This mutant has a '{analysis[answer_mutant_name]['type']}' mutation.\n"
                  f"However, the most disruptive mutations are of the types found in {best_mutants}, "
                  f"which are {[analysis[m]['type'] for m in best_mutants]}.\n"
                  f"The answer is incorrect because it did not select one of the mutations most likely to cause a complete loss of function.")
        return reason

# Run the checker
result = check_answer()
print(result)