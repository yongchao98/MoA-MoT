def explain_latitude_voc_effect():
    """
    Explains the theoretical effect of latitude on plant VOC alpha and beta diversity
    based on the plant-parasite arms race hypothesis.
    """

    # --- Core Principle ---
    print("Ecological Principle: The intensity of biotic interactions, particularly from parasites and herbivores, is highest in the tropics (low latitudes) and decreases towards the poles (high latitudes).")
    print("-" * 20)

    # --- Alpha Diversity (within-plant) ---
    print("Analysis of Alpha (α) Diversity (VOCs within a single plant):")
    print("1. At low latitudes (e.g., equator), high parasite pressure forces individual plants to produce a more complex mixture of defensive VOCs.")
    print("2. At high latitudes (e.g., 60°N), lower parasite pressure means a less diverse chemical arsenal is sufficient.")
    print("3. Therefore, as latitude decreases, α-diversity increases.")
    print("Conclusion for α-diversity: The relationship is NEGATIVE.")
    print("-" * 20)

    # --- Beta Diversity (between-plants at a site) ---
    print("Analysis of Beta (β) Diversity (VOC variation among plants at a site):")
    print("1. At low latitudes, high parasite pressure favors high chemical variation among plants. This prevents specialized parasites from easily locating and consuming their host plants.")
    print("2. At high latitudes, weaker parasite pressure and stronger abiotic filtering lead to more uniform chemical strategies among plants.")
    print("3. Therefore, as latitude decreases, β-diversity increases.")
    print("Conclusion for β-diversity: The relationship is NEGATIVE.")
    print("-" * 20)

    # --- Final Summary ---
    print("Final Summary:")
    print("The predicted direction of effect of latitude on VOC diversity is:")
    print("α-diversity: negative")
    print("β-diversity: negative")

# Execute the explanation
explain_latitude_voc_effect()