def analyze_insect_wings():
    """
    Analyzes three insect forewings to determine the trophic level of the family each belongs to.
    """
    print("Analyzing the trophic levels based on wing morphology:\n")

    # Analysis for Wing A
    print("--- Wing A ---")
    print("Identification: This wing belongs to the order Hymenoptera. The venation, characterized by a prominent pterostigma and a specific arrangement of closed cells, is typical of the superfamily Ichneumonoidea (e.g., Ichneumonidae or Braconidae).")
    print("Venation (Comstock-Needham): The venation is reduced. Key features are the Costa (C), a fused Subcosta (Sc) and Radius (R), a prominent pterostigma on the costal margin, and a distinct pattern of cells formed by the Media (M) and Cubitus (Cu) veins.")
    print("Trophic Level: Members of Ichneumonoidea are almost universally parasitoids, whose larvae feed on and kill a host arthropod.")
    trophic_level_A = "Parasitoid"
    print(f"Conclusion: Wing A represents a {trophic_level_A}.\n")

    # Analysis for Wing B
    print("--- Wing B ---")
    print("Identification: This wing is from the order Diptera, characteristic of the family Syrphidae (hoverflies).")
    print("Venation (Comstock-Needham): Key features include a 'spurious vein' (vena spuria) located between the Radius (R) and Media (M) veins, which is not a true, articulated vein. The radial veins have a distinctive curvature near the wing margin.")
    print("Trophic Level: While adult hoverflies are pollinators, many syrphid larvae are significant predators, especially of aphids. As the question asks for the highest possible trophic level, we classify this family as predatory.")
    trophic_level_B = "Predator"
    print(f"Conclusion: Wing B represents a {trophic_level_B}.\n")

    # Analysis for Wing C
    print("--- Wing C ---")
    print("Identification: This is a hardened, leathery forewing known as a tegmen. This structure is characteristic of the order Orthoptera (grasshoppers, katydids, etc.).")
    print("Venation (Comstock-Needham): The venation is primarily for structural support rather than aerodynamics. It consists of strong, roughly parallel longitudinal veins (C, Sc, R, M, Cu, A) connected by numerous crossveins.")
    print("Trophic Level: The overwhelming majority of Orthopterans are plant-feeders.")
    trophic_level_C = "Herbivore"
    print(f"Conclusion: Wing C represents a {trophic_level_C}.\n")

    # Final Summary
    print("--- Summary ---")
    print(f"A: {trophic_level_A}")
    print(f"B: {trophic_level_B}")
    print(f"C: {trophic_level_C}")
    print("\nThis combination matches answer choice E.")

analyze_insect_wings()