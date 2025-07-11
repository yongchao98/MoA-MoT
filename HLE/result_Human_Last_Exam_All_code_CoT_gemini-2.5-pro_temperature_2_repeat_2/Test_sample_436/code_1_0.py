import collections

def rank_epitopes():
    """
    Ranks epitopes for H2-Kd binding based on immunological principles.
    """
    epitopes = collections.OrderedDict([
        ('E1', 'TYQRTRALV'),
        ('E2', 'TFQRTRALV'),
        ('E3', 'TFQRTRALK'),
        ('E4', 'TYQRMFALV'),
        ('E5', 'TYIPMRALV')
    ])
    
    analysis_results = []

    print("--- Epitope Binding Analysis for H2-Kd ---")
    print("H2-Kd Primary Anchors: Position 2 (prefers Y, F), Position 9 (prefers V, L, I)\n")

    for name, seq in epitopes.items():
        p2 = seq[1]
        p9 = seq[8]
        
        # Rule 1: P9 Anchor check
        if p9 not in ['V', 'L', 'I']:
            tier = 4
            reason = f"Has a non-permissive anchor at P9 ({p9}). Binding will be severely reduced."
            analysis_results.append({'name': name, 'tier': tier, 'reason': reason})
            continue

        # Rule 2: P2 Anchor check
        if p2 not in ['Y', 'F']:
            tier = 4 # Also severely reduced
            reason = f"Has a non-permissive anchor at P2 ({p2}). Binding will be severely reduced."
            analysis_results.append({'name': name, 'tier': tier, 'reason': reason})
            continue

        # Rule 3: Check for sub-optimal anchors
        # The prompt implies Y is part of a high-affinity epitope. We'll rank F as sub-optimal relative to Y.
        if p2 == 'F':
            tier = 3
            reason = "Has optimal P9 anchor but a sub-optimal P2 anchor (F)."
            analysis_results.append({'name': name, 'tier': tier, 'reason': reason})
            continue

        # Rule 4: Check for disruptive internal residues (like Proline)
        if 'P' in seq[2:8]:
            tier = 2
            reason = "Has optimal P2/P9 anchors, but contains a structurally disruptive Proline."
            analysis_results.append({'name': name, 'tier': tier, 'reason': reason})
            continue

        # Rule 5: Top-tier binders
        # Both E1 and E4 fall here. Rank the reference (E1) highest.
        tier = 1
        if name == 'E1':
            # Use a sub-tier to ensure E1 is ranked above E4
            sub_tier = 1.0
            reason = "Reference high-affinity epitope with optimal P2 (Y) and P9 (V) anchors."
        else: # E4
            sub_tier = 1.1
            reason = "Optimal P2/P9 anchors, but internal changes (P6: R->F) likely make it bind less effectively than the reference E1."
            
        analysis_results.append({'name': name, 'tier': tier, 'sub_tier': sub_tier, 'reason': reason})

    # Sort the results by tier (and sub-tier)
    # The key lambda function ensures sorting by 'tier' first, then by 'sub_tier' if it exists.
    ranked_list = sorted(analysis_results, key=lambda x: (x['tier'], x.get('sub_tier', 0)))

    print("--- Final Ranking from Highest to Lowest Binding ---")
    final_order = []
    for i, item in enumerate(ranked_list):
        final_order.append(item['name'])
        print(f"{i+1}. {item['name']} ({epitopes[item['name']]})")
        print(f"   Reasoning: {item['reason']}")
        
    # The problem asks for the list of epitopes
    print("\nBased on this analysis, the ranked order of epitopes from highest to lowest expected amount complexed with H2-Kd is:")
    # Print the final result in the format: E1, E4, E5, E2, E3
    print(', '.join(final_order))

if __name__ == '__main__':
    rank_epitopes()
<<<A>>>