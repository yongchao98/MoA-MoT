import math

def analyze_proteomics_data():
    """
    Analyzes the volcano plot data to identify the protein with a reduced degradation rate.
    """
    # Approximate coordinates (LOG2(FOLD), -LOG10(P-VALUE)) from the plot
    proteins = {
        'A': (-4.8, 5.0),
        'B': (-3.5, 5.8),
        'C': (0.8, 3.0),
        'D': (3.2, 6.2),
        'E': (4.8, 3.5)
    }

    print("Step 1: Understanding the biological context.")
    print("A 'reduction in degradation rate' means the protein is broken down more slowly.")
    print("This leads to an accumulation of the protein, so its overall level should INCREASE.\n")

    print("Step 2: Interpreting the volcano plot.")
    print("We need to find a protein with:")
    print("  - A positive LOG2(FOLD) value (increase in protein level).")
    print("  - A high -LOG10(P-VALUE) (the increase is statistically significant).\n")

    print("Step 3: Evaluating each candidate protein.")
    
    best_candidate = None
    max_significance_for_increase = -1

    for protein, (log2_fold, neg_log_p_value) in proteins.items():
        fold_change = 2**log2_fold
        p_value = 10**(-neg_log_p_value)
        
        print(f"--- Protein {protein} ---")
        print(f"LOG2(FOLD) = {log2_fold:.2f}, -LOG10(P-VALUE) = {neg_log_p_value:.2f}")
        
        if log2_fold > 0:
            print(f"Result: This protein has INCREASED by a fold of {fold_change:.2f}.")
            print(f"The change is significant (p-value = {p_value:.2e}).")
            # Check if this is the best candidate so far among the increased proteins
            if neg_log_p_value > max_significance_for_increase:
                max_significance_for_increase = neg_log_p_value
                best_candidate = protein
        else:
            print(f"Result: This protein has DECREASED. This does not match the condition.\n")
            
    print("\nStep 4: Final Conclusion.")
    print(f"Proteins C, D, and E all show an increase.")
    print(f"Among these, Protein {best_candidate} shows a large increase combined with the highest statistical significance (-LOG10(P-VALUE) = {max_significance_for_increase:.2f}).")
    print(f"Therefore, Protein {best_candidate} is the most likely protein to have a reduced degradation rate.")

analyze_proteomics_data()