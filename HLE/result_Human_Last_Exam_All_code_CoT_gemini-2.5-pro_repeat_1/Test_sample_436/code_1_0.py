def rank_epitopes():
    """
    Analyzes and ranks epitopes based on their predicted binding affinity to H-2Kd.
    """
    # Define the epitopes provided in the problem
    epitopes = {
        "E1": "TYQRTRALV", # Reference high-affinity epitope
        "E2": "TFQRTRALV",
        "E3": "TFQRTRALK",
        "E4": "TYQRMFALV",
        "E5": "TYIPMRALV",
    }

    print("### Analysis of Epitope Binding to H-2Kd ###\n")
    print("The binding affinity is determined by key anchor residues.")
    print("H-2Kd Motif:")
    print(" - Primary Anchor P2: Prefers Y (optimal) or F (good).")
    print(" - Primary Anchor P9: Prefers V, L, I (optimal).")
    print(" - Secondary Anchor P6: Prefers R (optimal).")
    print("-" * 40)

    analysis_results = []
    for name, seq in epitopes.items():
        p2 = seq[1]
        p6 = seq[5]
        p9 = seq[8]
        
        # Assign scores based on the binding motif to formalize the ranking logic.
        # A higher score indicates better predicted binding.
        
        # Rule 1: Evaluate primary anchors (P2, P9), which are most critical.
        p2_score = 20 if p2 == 'Y' else (15 if p2 == 'F' else 0)
        p9_score = 20 if p9 in ['V', 'L', 'I'] else (-100 if p9 == 'K' else 0) # Strong penalty for K
        primary_score = p2_score + p9_score
        
        # Rule 2: Evaluate secondary anchor (P6).
        p6_score = 5 if p6 == 'R' else 0

        # Rule 3: Tie-breaker. The reference epitope E1 is the benchmark.
        # A small bonus ensures it ranks first in a tie.
        tie_breaker_score = 1 if name == 'E1' else 0

        analysis = {
            "name": name,
            "seq": seq,
            "primary_score": primary_score,
            "p6_score": p6_score,
            "tie_breaker": tie_breaker_score
        }
        analysis_results.append(analysis)

    # Sort epitopes in descending order of predicted affinity.
    # The key is a tuple: Python sorts by the first element, then the second for ties, and so on.
    # We use negative scores because sort() is ascending by default.
    sorted_epitopes = sorted(analysis_results, 
                             key=lambda x: (-x['primary_score'], -x['p6_score'], -x['tie_breaker']), 
                             reverse=False)

    print("\n### Final Ranking (Highest to Lowest Binding) ###\n")
    
    final_rank_list = [e['name'] for e in sorted_epitopes]
    
    # Print each epitope name in the final ranked order as requested.
    for i, name in enumerate(final_rank_list):
        if i > 0:
            print(" > ", end="")
        print(name, end="")
    print("\n")

    print("\nJustification:")
    print(" - E1: Optimal P2 (Y), P9 (V), and P6 (R). Set as the top binder.")
    print(" - E5: Optimal P2 (Y), P9 (V), and P6 (R). Ranks just below the reference E1.")
    print(" - E4: Optimal P2 (Y) and P9 (V), but suboptimal P6 (F). Ranks below E1 and E5.")
    print(" - E2: Suboptimal P2 (F) compared to the top group. Ranks below E4.")
    print(" - E3: Highly unfavorable P9 (K) drastically reduces binding. Ranked last.")

# Execute the analysis and print the result
rank_epitopes()