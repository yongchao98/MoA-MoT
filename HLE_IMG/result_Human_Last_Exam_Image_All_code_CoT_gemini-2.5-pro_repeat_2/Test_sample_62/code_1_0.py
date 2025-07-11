import math

def analyze_proteomics_data():
    """
    Analyzes volcano plot data to identify the protein with reduced degradation.
    """

    # Step 1: Approximate data points for the labeled proteins from the volcano plot.
    # The data is stored as a dictionary: {Protein: (LOG2(FOLD), -LOG10(P-VALUE))}
    proteins = {
        'A': (-4.8, 5.0),
        'B': (-3.8, 5.8),
        'C': (0.8, 3.1),
        'D': (3.2, 6.2),
        'E': (4.8, 3.5)
    }

    # Step 2: State the biological principle.
    print("Biological Rationale:")
    print("The level of a protein is a balance between its synthesis and degradation.")
    print("A 'reduction in its degradation rate' means the protein is destroyed more slowly.")
    print("This leads to an accumulation of the protein, so we expect its measured level to INCREASE.\n")

    print("Analysis of the Volcano Plot:")
    print("An INCREASE in protein level corresponds to a positive LOG2(FOLD) value (the x-axis).")
    print("An 'important' change corresponds to a high statistical significance, which is a high -LOG10(P-VALUE) (the y-axis).\n")

    # Step 3: Analyze each protein based on the principle.
    print("Evaluating each candidate protein:")
    candidates = {}
    for name, (log2_fold, neg_log10_p_value) in proteins.items():
        fold_change = 2**log2_fold if log2_fold >= 0 else -(2**(-log2_fold))
        p_value = 10**(-neg_log10_p_value)
        
        print(f"--- Protein {name} ---")
        print(f"Coordinates (LOG2 FOLD, -LOG10 P-VALUE): ({log2_fold}, {neg_log10_p_value})")
        
        if log2_fold > 0:
            print(f"Result: Protein level INCREASED approximately {fold_change:.2f}-fold.")
            print("This is consistent with a reduced degradation rate. It is a potential candidate.\n")
            candidates[name] = proteins[name]
        else:
            print(f"Result: Protein level DECREASED approximately {abs(fold_change):.2f}-fold.")
            print("This is NOT consistent with a reduced degradation rate. It is eliminated as a candidate.\n")

    # Step 4: Determine the best candidate from the list of possibilities.
    print("--- Final Conclusion ---")
    print("The candidates showing an increase are:", ", ".join(candidates.keys()))
    
    # Find the candidate with the highest statistical significance.
    best_candidate = ''
    max_p_value_score = -1
    for name, (log2_fold, neg_log10_p_value) in candidates.items():
        if neg_log10_p_value > max_p_value_score:
            max_p_value_score = neg_log10_p_value
            best_candidate = name
            
    print(f"Among these candidates, Protein {best_candidate} shows both a strong increase (LOG2 FOLD = {candidates[best_candidate][0]})")
    print(f"and the highest statistical significance (-LOG10 P-VALUE = {candidates[best_candidate][1]}).")
    print("This makes it the most likely protein to have an important reduction in its degradation rate.")

# Execute the analysis
analyze_proteomics_data()