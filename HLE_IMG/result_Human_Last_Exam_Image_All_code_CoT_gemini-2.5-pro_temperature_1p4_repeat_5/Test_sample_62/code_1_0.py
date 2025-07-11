import math

def analyze_proteomics_data():
    """
    Analyzes the volcano plot to identify the protein with a reduced degradation rate,
    and prints the reasoning.
    """
    print("Step 1: Determine the expected effect on the plot.")
    print("A reduced protein degradation rate leads to protein accumulation, causing an increase in its measured abundance.")
    print("On the volcano plot, an increase in abundance corresponds to a positive LOG2(FOLD) value (right side of the plot).\n")

    print("Step 2: Examine the labeled proteins on the plot.")
    # Approximate coordinates read from the graph
    protein_data = {
        'A': {'log2_fold': -4.8},
        'B': {'log2_fold': -3.8},
        'C': {'log2_fold': 0.8},
        'D': {'log2_fold': 3.2, 'neg_log10_p': 6.2},
        'E': {'log2_fold': 4.8, 'neg_log10_p': 3.5}
    }

    print("Proteins A and B have negative LOG2(FOLD) values, so their levels decreased.")
    print("Proteins C, D, and E have positive LOG2(FOLD) values, so their levels increased. They are our candidates.\n")

    print("Step 3: Select the best candidate among C, D, and E.")
    print("The term 'important reduction' suggests a result that is both large and statistically significant.")
    
    # Analysis of Protein D
    d_log2_fold = protein_data['D']['log2_fold']
    d_fold_change = math.pow(2, d_log2_fold)
    d_significance = protein_data['D']['neg_log10_p']
    
    print(f"\n- Protein D shows a large increase in abundance:")
    print(f"  LOG2(FOLD) = {d_log2_fold}")
    print(f"  This is calculated as: FOLD_CHANGE = 2^{d_log2_fold} = {d_fold_change:.1f}")
    print(f"  The fold change is approximately {d_fold_change:.1f}-fold.\n")
    print(f"- Protein D also shows the highest statistical significance on the plot:")
    print(f"  -LOG10(P-VALUE) = {d_significance}\n")

    print("Conclusion: Protein D has both a large fold-change and the highest statistical significance.")
    print("This makes it the most robust finding and the best candidate for the protein with an 'important' reduction in its degradation rate.")

analyze_proteomics_data()