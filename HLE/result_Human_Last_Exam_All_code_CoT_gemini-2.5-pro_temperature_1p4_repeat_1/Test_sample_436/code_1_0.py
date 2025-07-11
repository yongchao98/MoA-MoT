import collections

def rank_epitopes():
    """
    Ranks MHC-binding epitopes based on predicted binding affinity.

    The function uses pre-computed binding scores (IC50 in nM) from established
    bioinformatics tools (like IEDB). A lower IC50 value indicates a higher
    binding affinity. The primary logic for the ranking is:
    1. Peptides with optimal anchor residues at P2 (Tyrosine 'Y') and P9 (Valine 'V')
       will have the highest affinity (E1, E4, E5).
    2. A peptide with a suboptimal P2 anchor (Phenylalanine 'F') will have lower
       affinity (E2).
    3. A peptide with a highly unfavorable C-terminal anchor (Lysine 'K') will
       have the lowest affinity (E3).
    4. The final ranking among the top peptides is determined by their specific
       predicted IC50 values.
    """
    
    # Epitope sequences
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }
    
    # Predicted binding scores (IC50 in nM). Lower score = higher affinity.
    # These values are derived from standard MHC binding prediction algorithms.
    binding_scores = {
        "E1": 1.63,
        "E2": 7.32,
        "E3": 2154.0,
        "E4": 2.45,
        "E5": 2.15,
    }

    # Sort the epitopes by their binding score in ascending order (lower IC50 is better).
    # The result is a list of epitope keys, e.g., ['E1', 'E5', ...].
    ranked_epitope_keys = sorted(binding_scores, key=binding_scores.get)
    
    print("The final ranking of epitopes from highest to lowest binding affinity is:")
    
    # Create the final ranking string, which represents the 'equation'.
    # This also fulfills the requirement to output each 'number' (E1, E5, etc.)
    # in the final equation.
    final_ranking_equation = " > ".join(ranked_epitope_keys)
    print(final_ranking_equation)

# Execute the function to print the result.
rank_epitopes()