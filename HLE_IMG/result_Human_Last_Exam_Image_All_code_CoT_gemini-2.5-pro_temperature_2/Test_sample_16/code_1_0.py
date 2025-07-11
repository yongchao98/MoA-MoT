def solve_insect_wings():
    """
    Analyzes the provided insect wings to determine the trophic level of each and select the correct answer.
    """
    
    # Define the descriptions based on entomological analysis
    descriptions = {
        "A": {
            "venation": "The forewing is membranous, with venation characteristic of the Order Hymenoptera. A key feature is the prominent dark pterostigma on the costal (leading) edge. The arrangement of closed cells formed by the Radius, Media, and Cubitus veins and their crossveins is typical for a parasitic wasp.",
            "identification": "The wing belongs to a parasitic wasp, likely from the superfamily Ichneumonoidea (e.g., family Ichneumonidae).",
            "trophic_level": "Parasitoid"
        },
        "B": {
            "venation": "The forewing is broad, membranous, and displays a complex, net-like (reticulate) venation characteristic of the Order Neuroptera. There are numerous costal crossveins and extensive branching of all major longitudinal veins (Radius, Media, Cubitus), creating a high density of cells across the entire wing.",
            "identification": "The wing belongs to a lacewing (e.g., family Chrysopidae or Hemerobiidae).",
            "trophic_level": "Predator"
        },
        "C": {
            "venation": "This is a hardened, leathery forewing known as a tegmen, characteristic of the Order Orthoptera. The longitudinal veins run largely parallel, and the areas between them are filled with a dense network of crossveins, giving it a reticulated appearance.",
            "identification": "The wing belongs to an Orthopteran, such as a grasshopper or katydid.",
            "trophic_level": "Herbivore"
        }
    }
    
    # Print the thorough description
    print("Detailed analysis of each wing:\n")
    for wing_id in ["A", "B", "C"]:
        print(f"--- Wing {wing_id} ---")
        print(f"Venation Description: {descriptions[wing_id]['venation']}")
        print(f"Identification and Trophic Level: {descriptions[wing_id]['identification']} Both larvae and adults of this group are primarily known for their feeding habit, which places them in the **{descriptions[wing_id]['trophic_level']}** trophic level.")
        print("")

    # Print the final conclusion
    print("\n--- Conclusion ---")
    print("Based on the wing morphology:")
    print(f"Wing A belongs to a {descriptions['A']['trophic_level']}.")
    print(f"Wing B belongs to a {descriptions['B']['trophic_level']}.")
    print(f"Wing C belongs to a {descriptions['C']['trophic_level']}.")
    print("\nThis combination matches choice E.")

    # Print the final answer in the required format
    print("<<<E>>>")

solve_insect_wings()