import math

def analyze_proteomics_data():
    """
    Analyzes protein data from a volcano plot to identify the protein
    with a reduced degradation rate.
    """

    # Approximate coordinates (LOG2(FOLD), -LOG10(P-VALUE)) from the plot
    proteins = {
        'A': {'log2_fold': -4.8, 'neg_log10_p_value': 5.0},
        'B': {'log2_fold': -3.6, 'neg_log10_p_value': 5.8},
        'C': {'log2_fold': 0.8, 'neg_log10_p_value': 3.2},
        'D': {'log2_fold': 3.2, 'neg_log10_p_value': 6.1},
        'E': {'log2_fold': 4.8, 'neg_log10_p_value': 3.5}
    }

    print("Step 1: Understanding the biological context.")
    print("A reduced protein degradation rate means the protein is destroyed more slowly.")
    print("This leads to an accumulation of the protein, meaning its abundance increases.\n")

    print("Step 2: Relating the context to the volcano plot.")
    print("An increase in protein abundance corresponds to a positive 'LOG2 (FOLD)' value (up-regulation).")
    print("An 'important' change implies high statistical significance, which corresponds to a high '-LOG10 (P-VALUE)'.\n")

    print("Step 3: Filtering for up-regulated proteins (LOG2_FOLD > 0).")
    upregulated_proteins = {}
    for name, data in proteins.items():
        if data['log2_fold'] > 0:
            upregulated_proteins[name] = data
            print(f"Protein {name} is up-regulated with LOG2(FOLD) = {data['log2_fold']:.2f}. It is a candidate.")
        else:
            print(f"Protein {name} is down-regulated with LOG2(FOLD) = {data['log2_fold']:.2f}. It is not a candidate.")
    print("\nCandidates are:", list(upregulated_proteins.keys()), "\n")

    print("Step 4: Finding the most statistically significant protein among the candidates.")
    most_significant_protein = None
    max_p_value = -1

    for name, data in upregulated_proteins.items():
        print(f"Checking candidate {name}: -LOG10(P-VALUE) = {data['neg_log10_p_value']:.2f}")
        if data['neg_log10_p_value'] > max_p_value:
            max_p_value = data['neg_log10_p_value']
            most_significant_protein = name

    print(f"\nThe protein with the highest statistical significance (-LOG10(P-VALUE) = {max_p_value:.2f}) is Protein {most_significant_protein}.")
    
    log2_fold_d = proteins['D']['log2_fold']
    fold_change = 2**log2_fold_d
    neg_log10_p_value_d = proteins['D']['neg_log10_p_value']
    p_value = 10**(-neg_log10_p_value_d)

    print("\n--- Conclusion ---")
    print(f"Protein D is the most likely candidate because it is both significantly up-regulated (its abundance increased ~{fold_change:.1f}-fold) and the change is the most statistically significant (p-value ≈ {p_value:.1e}).")
    print("This profile is consistent with a protein accumulating due to a reduced degradation rate.")
    print("\nThe question asks for the protein. The answer is Protein D.")
    print("The corresponding answer choice is A.")


analyze_proteomics_data()