import math

def analyze_proteomics_data():
    """
    Analyzes the volcano plot data to identify the protein with a reduced degradation rate.
    """
    # Estimated data points from the volcano plot for the labeled proteins
    proteins = {
        'A': {'log2_fold_change': -4.8, 'neg_log10_p_value': 5.0},
        'B': {'log2_fold_change': -3.6, 'neg_log10_p_value': 5.8},
        'C': {'log2_fold_change': 0.8, 'neg_log10_p_value': 3.2},
        'D': {'log2_fold_change': 3.2, 'neg_log10_p_value': 6.2},
        'E': {'log2_fold_change': 4.8, 'neg_log10_p_value': 3.5}
    }

    print("Step 1: Analyzing the biological context.")
    print("A reduced degradation rate leads to protein accumulation, which means an INCREASE in the protein's level.")
    print("-" * 30)

    print("Step 2: Interpreting the volcano plot.")
    print("We are looking for a protein with a significant INCREASE in level.")
    print("This means we need a point with:")
    print("  - A large POSITIVE LOG2(FOLD) value (x-axis).")
    print("  - A large -LOG10(P-VALUE) value (y-axis for significance).")
    print("-" * 30)

    print("Step 3: Evaluating each labeled protein.")
    best_candidate = None
    max_significance_in_upregulated = -1

    for name, data in proteins.items():
        log2_fold = data['log2_fold_change']
        p_value_metric = data['neg_log10_p_value']
        
        # Calculate the actual fold change from the log2 value
        fold_change = 2**log2_fold
        
        if log2_fold > 0:
            print(f"Protein {name}: LOG2(FOLD) = {log2_fold}, -LOG10(P-VALUE) = {p_value_metric}.")
            print(f"  - Interpretation: Level increased by a fold change of {fold_change:.2f}. This is a potential candidate.")
            if p_value_metric > max_significance_in_upregulated:
                max_significance_in_upregulated = p_value_metric
                best_candidate = name
        else:
            print(f"Protein {name}: LOG2(FOLD) = {log2_fold}, -LOG10(P-VALUE) = {p_value_metric}.")
            print(f"  - Interpretation: Level decreased. This is not a candidate.")
    print("-" * 30)
    
    print("Step 4: Final Conclusion.")
    print(f"Among the up-regulated proteins (C, D, E), Protein D has a large fold change (LOG2(FOLD) = {proteins['D']['log2_fold_change']}) and the highest statistical significance (-LOG10(P-VALUE) = {proteins['D']['neg_log10_p_value']}).")
    print("This makes it the most 'important' change consistent with a reduced degradation rate.")
    print(f"\nThe correct protein is {best_candidate}.")

# Run the analysis
analyze_proteomics_data()