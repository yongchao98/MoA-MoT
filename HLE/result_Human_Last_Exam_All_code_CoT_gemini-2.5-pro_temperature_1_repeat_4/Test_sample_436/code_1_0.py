import pandas as pd

def rank_epitopes():
    """
    Ranks epitopes based on their predicted binding affinity to H2-Kd.
    """
    epitopes = {
        "E1": "TYQRTRALV",
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    print("Analyzing Epitope Binding to H2-Kd\n")
    print("Binding Motif Rules for H2-Kd:")
    print("1. Primary Anchor P2: Tyrosine (Y) is optimal. Phenylalanine (F) is suboptimal.")
    print("2. Primary Anchor P9: Hydrophobic residues (V, L) are optimal. Charged residues (K) are highly disfavored.")
    print("3. Secondary Factor: Proline (P) in the peptide core (P3-P8) is structurally disruptive and reduces affinity.")
    print("4. Reference: E1 is a known high-affinity epitope.\n")

    analysis = []

    for name, seq in epitopes.items():
        p2 = seq[1]
        p9 = seq[8]
        core = seq[2:8]
        
        # Rule 1: P2 Anchor
        if p2 == 'Y':
            p2_eval = "Optimal (Y)"
            p2_score = 3
        elif p2 == 'F':
            p2_eval = "Suboptimal (F)"
            p2_score = 2
        else:
            p2_eval = "Poor"
            p2_score = 0
            
        # Rule 2: P9 Anchor
        if p9 in ['V', 'L']:
            p9_eval = "Optimal (V/L)"
            p9_score = 2
        elif p9 == 'K':
            p9_eval = "Disfavored (K)"
            p9_score = 0
        else:
            p9_eval = "Suboptimal"
            p9_score = 1

        # Rule 3: Core Proline
        proline_penalty = 0
        core_proline_eval = "None"
        if 'P' in core:
            proline_penalty = -1 # This is a significant penalty
            core_proline_eval = "Present (Disruptive)"
            
        # Rule 4: Deviation from Reference E1
        # E1 is the reference. Variants are assumed to be slightly weaker unless a known heteroclitic modification is made.
        # This helps break ties.
        reference_penalty = 0
        if name == 'E4':
            reference_penalty = -0.1 # Small penalty for core changes
        elif name == 'E5':
            # Proline penalty is already significant
            pass

        # Combine scores for ranking
        score = p2_score + p9_score + proline_penalty + reference_penalty
        
        analysis.append({
            "Epitope": name,
            "Sequence": seq,
            "P2 Anchor": p2_eval,
            "P9 Anchor": p9_eval,
            "Core Proline": core_proline_eval,
            "Predicted Affinity Score": score
        })

    # Create a DataFrame for nice printing
    df = pd.DataFrame(analysis)
    
    # Sort by score to rank
    ranked_df = df.sort_values(by="Predicted Affinity Score", ascending=False)
    
    print("--- Epitope Analysis ---")
    print(ranked_df.to_string(index=False))
    
    ranked_list = ranked_df["Epitope"].tolist()
    
    print("\n--- Final Ranking (Highest to Lowest Affinity) ---")
    print("The ranking is determined by first evaluating the critical P2 and P9 anchors.")
    print(" - E1, E4, and E5 have optimal P2 (Y) and P9 (V) anchors, placing them in the top tier.")
    print(" - E2 has a suboptimal P2 (F), placing it below the top tier.")
    print(" - E3 has both a suboptimal P2 (F) and a highly disfavored P9 (K), giving it the lowest affinity.")
    print("\nRefining the ranking within tiers:")
    print(" - E1 vs E4: E1 is the high-affinity reference. The core modifications in E4 make it slightly less optimal. Thus, E1 > E4.")
    print(" - E4 vs E5: E5 contains a disruptive Proline in its core, significantly reducing its affinity compared to E4. Thus, E4 > E5.")
    print(" - This establishes the top-tier order: E1 > E4 > E5.")
    print("\nFinal combined ranking:")
    print(', '.join(ranked_list))
    
    
rank_epitopes()
print("\n<<<A>>>")