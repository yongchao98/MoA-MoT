import math

def analyze_protein_dynamics():
    """
    Analyzes protein data from a volcano plot to identify a protein
    with a reduced degradation rate.
    """

    # Approximate coordinates of the labeled proteins from the volcano plot.
    # 'log2_fold': LOG2(FOLD CHANGE)
    # 'p_value': The corresponding p-value calculated from -LOG10(P-VALUE)
    proteins = {
        'A': {'log2_fold': -4.8, 'neg_log10_p': 5.0},
        'B': {'log2_fold': -3.6, 'neg_log10_p': 5.8},
        'C': {'log2_fold': 0.8, 'neg_log10_p': 3.2},
        'D': {'log2_fold': 3.2, 'neg_log10_p': 6.2},
        'E': {'log2_fold': 4.8, 'neg_log10_p': 3.6},
    }

    print("### Analysis of Protein Regulation ###\n")
    print("A reduction in protein degradation leads to protein accumulation, resulting in up-regulation.")
    print("Up-regulated proteins have a positive LOG2(FOLD) value.\n")

    candidates = []
    for protein, data in proteins.items():
        log2_fold = data['log2_fold']
        neg_log10_p = data['neg_log10_p']
        fold_change = 2**log2_fold
        
        if log2_fold > 0:
            regulation = "Up-regulated"
            conclusion = "This is a possible candidate."
            candidates.append(protein)
        else:
            regulation = "Down-regulated"
            conclusion = "This is not a candidate."
        
        print(f"--- Protein {protein} ---")
        print(f"LOG2(FOLD) = {log2_fold:.2f}, -LOG10(P-VALUE) = {neg_log10_p:.2f}")
        print(f"This protein is {regulation}.")
        print(f"Conclusion: {conclusion}\n")

    print("### Final Determination ###")
    print(f"The candidates are: {', '.join(candidates)}.")
    print("We are looking for a protein with an 'important' change. This implies both a large fold change and high statistical significance (high -LOG10(P-VALUE)).")

    # Find the candidate with the highest statistical significance.
    best_candidate = ''
    max_significance = -1
    for p_name in candidates:
      if proteins[p_name]['neg_log10_p'] > max_significance:
        max_significance = proteins[p_name]['neg_log10_p']
        best_candidate = p_name

    print(f"Comparing the candidates ({', '.join(candidates)}), Protein {best_candidate} has the highest statistical significance (-LOG10(P-VALUE) = {max_significance:.2f}) and a large fold change.")
    print("Therefore, Protein D is the most likely protein to have a reduced degradation rate.")

analyze_protein_dynamics()