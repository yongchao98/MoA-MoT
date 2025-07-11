def analyze_wings():
    """
    Analyzes the provided insect wings to determine the family and trophic level for each.
    """
    print("Analyzing Wing A:")
    print("  - Identification: This wing belongs to the order Hymenoptera, likely from the superfamily Ichneumonoidea (e.g., family Braconidae or Ichneumonidae).")
    print("  - Venation Description (Comstock-Needham System): The wing shows a prominent pterostigma on the costal margin. The venation is somewhat reduced, with a closed marginal (radial) cell and several submarginal and discal cells in the center of the wing.")
    print("  - Trophic Level: The families in Ichneumonoidea are overwhelmingly parasitoids, laying their eggs on or in other insects which serve as hosts for the developing larvae. Therefore, the highest trophic level is Parasitoid.\n")

    print("Analyzing Wing B:")
    print("  - Identification: This wing belongs to the order Diptera, family Syrphidae (hoverflies).")
    print("  - Venation Description (Comstock-Needham System): A key diagnostic feature is the vena spuria ('spurious vein'), a non-articulated, vein-like thickening between the Radius (R) and Media (M) veins. The R4+5 vein curves anteriorly, and the anal cell is distinctly long and closed before the wing margin.")
    print("  - Trophic Level: While adults are pollinators, the larvae of many hoverfly species are active predators of soft-bodied insects like aphids. The highest trophic level for this family is Predator.\n")

    print("Analyzing Wing C:")
    print("  - Identification: This wing is a hemelytron from the order Hemiptera, suborder Heteroptera (true bugs), likely from the superfamily Lygaeoidea (seed bugs).")
    print("  - Venation Description (Comstock-Needham System): The wing is divided into a coriaceous (leathery) base and a membranous apex. The venation is primarily visible in the membrane, which shows a pattern of several longitudinal veins forming a few closed cells.")
    print("  - Trophic Level: The majority of species in this group are phytophagous, feeding on plant tissues or seeds. Therefore, the trophic level is Herbivore.\n")

    print("Summary of Trophic Levels:")
    print("  - A: Parasitoid")
    print("  - B: Predator")
    print("  - C: Herbivore")

    print("\nThis corresponds to answer choice E.")

if __name__ == '__main__':
    analyze_wings()
    print("<<<E>>>")