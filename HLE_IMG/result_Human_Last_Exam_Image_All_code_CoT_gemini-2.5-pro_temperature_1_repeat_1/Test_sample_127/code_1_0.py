def solve_chemistry_rotation():
    """
    Analyzes the rotation of phenyl rings in the given iridium complex.
    """
    print("Analyzing the free rotation of phenyl rings G, H, J, K, and L.")
    print("-" * 30)

    # Dictionary to store the analysis for each ring
    analysis = {
        'H': "Attached at position 4 of ring A (para to the Ir-C bond). It is far from the metal center and points away from the complex. It can rotate freely.",
        'J': "Attached at position 4 of ring C (para to the Ir-C bond). Similar to ring H, it is far from the metal center and can rotate freely.",
        'G': "Attached at position 6 of pyridine ring D. It is on the periphery of the complex and is not sterically locked by the metal center. It can rotate freely.",
        'K': "Attached at position 6 of ring A (ortho to the Ir-C bond). It is very close to the Iridium center, causing severe steric hindrance that prevents free rotation.",
        'L': "Attached at position 6 of ring C (ortho to the Ir-C bond). Similar to ring K, its proximity to the Iridium center prevents free rotation."
    }

    freely_rotating_rings = []
    restricted_rings = []

    for ring, reason in analysis.items():
        print(f"Ring {ring}: {reason}")
        if "can rotate freely" in reason:
            freely_rotating_rings.append(ring)
        else:
            restricted_rings.append(ring)
    
    freely_rotating_rings.sort()

    print("-" * 30)
    print("Conclusion:")
    print(f"The rings that can rotate freely are: {', '.join(freely_rotating_rings)}")
    print(f"The rings with restricted rotation are: {', '.join(sorted(restricted_rings))}")
    print("\nMatching this with the answer choices, the correct option includes G, J, and H.")

solve_chemistry_rotation()
<<<E>>>