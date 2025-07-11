def analyze_persistence_of_fragments():
    """
    Analyzes the interplay between natural selection and genetic drift
    in determining the fate of small, non-functional genomic fragments.
    """
    # The selection coefficient (s) represents the fitness advantage of
    # deleting a small, useless genomic fragment.
    # It's a small positive value because replicating DNA costs energy.
    s = 0.0001
    print(f"Let's assume the selective advantage 's' of deleting a fragment is: {s}\n")

    # We will test a range of effective population sizes (Ne).
    # Endosymbionts and parasites often have small Ne.
    population_sizes = [100, 1000, 10000, 1000000]

    print("Analyzing the dominant evolutionary force for different population sizes:")
    print("-" * 65)

    for Ne in population_sizes:
        # The strength of genetic drift is inversely proportional to the
        # effective population size (Ne). We can represent its magnitude
        # relative to selection as 1/Ne.
        drift_strength_threshold = 1.0 / Ne

        print(f"Population Size (Ne): {Ne}")
        print(f"  - Selective Advantage (s) of deletion: {s}")
        print(f"  - Genetic Drift Strength Threshold (1/Ne): {drift_strength_threshold:.6f}")

        # If s is not much larger than 1/Ne, drift will be the dominant force.
        if s < drift_strength_threshold:
            dominance = "DRIFT"
            consequence = "Fragment is likely to persist due to random chance. Selection is too weak to be effective."
        else:
            dominance = "SELECTION"
            consequence = "Fragment is likely to be efficiently removed by natural selection."

        print(f"  - Dominant Force: {dominance}")
        print(f"  - Consequence: {consequence}\n")

    print("Conclusion: In small populations, genetic drift is strong and overpowers weak selection,")
    print("leading to the persistence of slightly deleterious genomic fragments. Therefore, the strength")
    print("of genetic drift is the primary factor influencing their persistence.")

# Run the analysis
analyze_persistence_of_fragments()