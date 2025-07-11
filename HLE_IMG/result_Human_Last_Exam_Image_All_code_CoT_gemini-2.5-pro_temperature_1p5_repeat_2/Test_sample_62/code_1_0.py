import math

def solve_proteomics_puzzle():
    """
    This script analyzes the data from the volcano plot to identify the correct protein.
    """
    # Approximate coordinates from the volcano plot
    proteins = {
        'A': {'log2_fold': -4.7, 'neg_log10_p_value': 4.9},
        'B': {'log2_fold': -3.7, 'neg_log10_p_value': 5.8},
        'C': {'log2_fold': 0.8, 'neg_log10_p_value': 3.0},
        'D': {'log2_fold': 3.2, 'neg_log10_p_value': 6.2},
        'E': {'log2_fold': 4.8, 'neg_log10_p_value': 3.5}
    }

    print("Step 1: Understand the biological context.")
    print("A 'reduction in degradation rate' means the protein breaks down slower, leading to its accumulation.")
    print("This means we are looking for an up-regulated protein.")
    print("-" * 30)

    print("Step 2: Relate the context to the plot axes.")
    print("Up-regulation corresponds to a positive LOG2 (FOLD) value (the x-axis).")
    print("An 'important' change corresponds to high statistical significance, which is a high -LOG10 (P-VALUE) (the y-axis).")
    print("-" * 30)
    
    print("Step 3: Evaluate each candidate based on the LOG2 (FOLD) value.")
    upregulated_proteins = {}
    for name, values in proteins.items():
        if values['log2_fold'] > 0:
            print(f"Protein {name} is UP-regulated (LOG2 Fold = {values['log2_fold']}). It is a potential candidate.")
            upregulated_proteins[name] = values
        else:
            print(f"Protein {name} is DOWN-regulated (LOG2 Fold = {values['log2_fold']}). It is eliminated.")
    print("-" * 30)

    print("Step 4: Find the most 'important' candidate among the up-regulated proteins.")
    print("We will find the protein with the highest statistical significance (highest -LOG10(P-VALUE)).")
    
    best_candidate = ''
    max_p_value = -1

    for name, values in upregulated_proteins.items():
        if values['neg_log10_p_value'] > max_p_value:
            max_p_value = values['neg_log10_p_value']
            best_candidate = name

    print(f"\nThe candidates are {list(upregulated_proteins.keys())}.")
    print(f"Protein {best_candidate} has the highest -LOG10(P-VALUE) of {max_p_value}.")
    fold_change = 2**proteins[best_candidate]['log2_fold']
    p_value = 10**(-proteins[best_candidate]['neg_log10_p_value'])
    print(f"This corresponds to an abundance increase of approximately {fold_change:.1f}-fold and a p-value of {p_value:.1e}.")
    print("\nConclusion: Protein D is the most significantly up-regulated protein, making it the best match for the description.")

solve_proteomics_puzzle()