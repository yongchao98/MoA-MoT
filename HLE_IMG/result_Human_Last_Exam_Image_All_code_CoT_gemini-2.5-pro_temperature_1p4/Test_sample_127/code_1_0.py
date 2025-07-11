def analyze_ring_rotation():
    """
    Analyzes the free rotation of phenyl rings in the given Iridium complex
    based on steric hindrance.
    """

    # Rings are classified based on their attachment point.
    # Rings attached to a cyclometalated carbon are considered restricted.
    # Other peripheral rings are considered free to rotate.

    all_rings = ['G', 'H', 'J', 'K', 'L']
    
    # Rings attached to carbons directly bonded to the Iridium center (metallated carbons).
    # These carbons are C-6 of ring A and C-6 of ring C.
    rings_on_metallated_carbons = ['K', 'L']
    
    freely_rotating_rings = []

    print("Analysis of Phenyl Ring Rotation:")
    print("=" * 35)

    for ring in sorted(all_rings):
        if ring in rings_on_metallated_carbons:
            print(f"Ring {ring}: Attached to a carbon directly bonded to the Ir metal center.")
            print("           This creates severe steric hindrance, restricting rotation.")
        else:
            print(f"Ring {ring}: Attached to the ligand periphery, away from the congested metal center.")
            print("           The single C-C bond allows for free rotation.")
            freely_rotating_rings.append(ring)
        print("-" * 35)

    print("\nConclusion:")
    print("The rings that can rotate freely are those not attached to the highly congested metalation sites.")
    
    # Sort for consistent output
    freely_rotating_rings.sort()
    
    # Use print() to display the final result clearly
    print("\nFreely rotating rings are:")
    for ring in freely_rotating_rings:
        print(f"- {ring}")

    # Combine the letters for the final answer
    final_answer_set = ", ".join(freely_rotating_rings)
    print(f"\nFinal Answer Set: {{{final_answer_set}}}")

# Run the analysis
analyze_ring_rotation()
