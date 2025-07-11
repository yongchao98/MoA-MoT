import math

def analyze_proteomics_data():
    """
    Analyzes proteomics data from a volcano plot to identify the protein
    with a reduced degradation rate.
    """
    # Estimated data points from the volcano plot
    # Format: {'Protein_Name': {'log2_fold': value, 'neg_log10_p': value}}
    protein_data = {
        'A': {'log2_fold': -4.8, 'neg_log10_p': 5.0},
        'B': {'log2_fold': -3.6, 'neg_log10_p': 5.8},
        'C': {'log2_fold': 0.8, 'neg_log10_p': 3.0},
        'D': {'log2_fold': 3.2, 'neg_log10_p': 6.2},
        'E': {'log2_fold': 4.8, 'neg_log10_p': 3.4}
    }

    print("Step 1: Analyzing the biological problem.")
    print("A reduced degradation rate leads to protein accumulation.")
    print("This means we are looking for a protein with a positive Log2(Fold Change), indicating it is up-regulated.\n")

    print("Step 2: Filtering proteins based on Log2(Fold Change).")
    up_regulated_proteins = {}
    for protein, data in protein_data.items():
        if data['log2_fold'] > 0:
            up_regulated_proteins[protein] = data
            print(f"- Protein {protein} is a candidate. (Log2 Fold = {data['log2_fold']})")
        else:
            print(f"- Protein {protein} is NOT a candidate. (Log2 Fold = {data['log2_fold']})")

    print("\nStep 3: Identifying the most 'important' change among candidates.")
    print("An 'important' change implies both a significant fold change and high statistical confidence (high -Log10(P-value)).\n")

    best_candidate = ''
    max_significance = -1

    for protein, data in up_regulated_proteins.items():
        log2_fold = data['log2_fold']
        neg_log10_p = data['neg_log10_p']
        
        # Determine the protein with the highest significance among the up-regulated ones
        if neg_log10_p > max_significance:
            max_significance = neg_log10_p
            best_candidate = protein
            
    print("Step 4: Final Conclusion.")
    # Calculate the values for the best candidate to show the final equation
    best_data = protein_data[best_candidate]
    fold_change = math.pow(2, best_data['log2_fold'])
    p_value = math.pow(10, -best_data['neg_log10_p'])
    
    print(f"Among the candidates (C, D, E), Protein D shows the highest statistical significance.")
    print(f"For Protein D: -Log10(P-value) = {best_data['neg_log10_p']:.1f} and Log2(Fold) = {best_data['log2_fold']:.1f}.")
    print(f"This translates to a p-value of {p_value:.2e} and a fold-change of {fold_change:.1f}-fold.")
    print("This combination of high significance and strong up-regulation makes it the best answer.")
    
    print(f"\nThe identified protein is: {best_candidate}")

# Execute the analysis
analyze_proteomics_data()