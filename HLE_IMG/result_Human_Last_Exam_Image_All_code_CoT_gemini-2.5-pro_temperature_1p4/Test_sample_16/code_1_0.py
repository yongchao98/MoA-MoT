def analyze_wings():
    """
    Analyzes three insect wings to determine the trophic level of each associated insect.
    """
    print("Analysis of Insect Forewings:\n")

    # Wing A Analysis
    print("--- Wing A ---")
    print("Description: This is a membranous wing with a prominent pterostigma (darkened spot on the leading edge).")
    print("Venation (Comstock-Needham System):")
    print(" - The venation is well-developed with a system of closed cells.")
    print(" - The Radius (R) vein encloses a distinct radial cell.")
    print(" - The Media (M) and Cubitus (Cu) veins are branched and connected by crossveins to form multiple cubital and discal cells.")
    print("Identification: This venation pattern is characteristic of a parasitic wasp (Order: Hymenoptera, likely Superfamily: Ichneumonoidea).")
    print("Trophic Level: Parasitoid\n")

    # Wing B Analysis
    print("--- Wing B ---")
    print("Description: This is a broad, membranous wing with a dense, complex network of veins.")
    print("Venation (Comstock-Needham System):")
    print(" - A key feature is the numerous simple crossveins along the leading (costal) margin.")
    print(" - The main longitudinal veins, especially the Radius (R) and Media (M), show extensive branching towards the wing tip (apical twigging).")
    print("Identification: This type of venation is classic for the Order Neuroptera (e.g., lacewings, antlions).")
    print("Trophic Level: Predator\n")

    # Wing C Analysis
    print("--- Wing C ---")
    print("Description: This is a leathery, hardened forewing known as a tegmen.")
    print("Venation (Comstock-Needham System):")
    print(" - The veins (R, M, Cu, etc.) are embedded within the leathery membrane and form a dense, net-like (reticulate) pattern for structural support.")
    print(" - This structure protects the delicate hind wing at rest.")
    print("Identification: Tegmina of this shape are characteristic of the Order Orthoptera (grasshoppers and katydids).")
    print("Trophic Level: Herbivore\n")

    # Final Conclusion
    print("--- Summary ---")
    print("A: Parasitoid")
    print("B: Predator")
    print("C: Herbivore")
    print("\nThis corresponds to Answer Choice E.")

analyze_wings()