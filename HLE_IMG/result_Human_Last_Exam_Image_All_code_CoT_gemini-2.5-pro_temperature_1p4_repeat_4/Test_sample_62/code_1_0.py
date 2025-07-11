import math

def analyze_proteomics_data():
    """
    Analyzes protein data from a volcano plot to identify the protein with a reduced degradation rate.
    """
    # Approximate coordinates of the labeled proteins from the volcano plot
    proteins = {
        'A': {'log2_fold': -4.8, 'neg_log10_p_value': 5.0},
        'B': {'log2_fold': -3.8, 'neg_log10_p_value': 5.8},
        'C': {'log2_fold': 0.8, 'neg_log10_p_value': 3.1},
        'D': {'log2_fold': 3.2, 'neg_log10_p_value': 6.2},
        'E': {'log2_fold': 4.8, 'neg_log10_p_value': 3.5}
    }

    print("Step 1: Understanding the Biological Problem")
    print("A reduction in a protein's degradation rate causes it to accumulate.")
    print("This means we are looking for a protein that is significantly *upregulated*.\n")

    print("Step 2: Interpreting the Volcano Plot")
    print("Upregulated proteins have a positive 'LOG2 (FOLD)' value (right side of the plot).")
    print("Statistically significant changes have a high '-LOG10 (P-VALUE)' (top of the plot).\n")

    print("Step 3: Analyzing Each Labeled Protein")
    upregulated_proteins = []
    for name, data in proteins.items():
        log2_fold = data['log2_fold']
        p_val_log = data['neg_log10_p_value']
        
        status = "Upregulated" if log2_fold > 0 else "Downregulated"
        print(f"Protein {name}: LOG2(FOLD) = {log2_fold}, -LOG10(P-VALUE) = {p_val_log}. Status: {status}")

        if log2_fold > 0:
            upregulated_proteins.append((name, data))
    
    print("\nStep 4: Identifying the Best Candidate")
    if not upregulated_proteins:
        print("No upregulated proteins found among the labeled points.")
        return

    # Find the upregulated protein with the highest significance
    # The highest significance corresponds to the highest -log10(p-value)
    best_candidate = max(upregulated_proteins, key=lambda p: p[1]['neg_log10_p_value'])
    
    candidate_name = best_candidate[0]
    candidate_data = best_candidate[1]
    
    print("The candidate proteins are those that are upregulated (C, D, E).")
    print("Among these, we look for the one with the highest statistical significance (highest y-value).")
    print(f"Protein {candidate_name} is the most significantly upregulated protein with a -LOG10(P-VALUE) of {candidate_data['neg_log10_p_value']}.")
    print("\nConclusion: Protein D best fits the description of a protein with an important reduction in its degradation rate.")


analyze_proteomics_data()