import textwrap

def analyze_insect_wings():
    """
    Analyzes three insect wings to determine the trophic level of the family each belongs to.
    """

    # --- Wing A Analysis ---
    wing_a_description = {
        "Identification": "Family Ichneumonidae (Ichneumon Wasps), Order Hymenoptera",
        "Venation (Comstock-Needham System)": [
            "Pterostigma: Prominent and heavily sclerotized, located on the costal (leading) margin.",
            "Veins C, Sc, R: Costa (C), Subcosta (Sc), and Radius (R) are present at the anterior of the wing.",
            "Submarginal Cells: There are typically two or three submarginal cells. A key feature is the second submarginal cell (also called the areolet), which is small and often pentagonal, resembling a 'horse head'. This feature is highly characteristic of Ichneumonidae.",
            "Other Cells: Numerous closed cells are formed by the Media (M) and Cubitus (Cu) veins, indicating a relatively complete venation pattern for a wasp."
        ],
        "Trophic Level": "Parasitoid. Ichneumonid wasps are primary parasitoids, with their larvae developing on or inside other arthropods (mostly other insect larvae), ultimately killing their host."
    }

    # --- Wing B Analysis ---
    wing_b_description = {
        "Identification": "Family Asilidae (Robber Flies), Order Diptera",
        "Venation (Comstock-Needham System)": [
            "Overall Shape: Typical fly wing, broader than the wasp wing in A.",
            "Veins: The venation is characteristic of the Brachycera suborder of flies.",
            "Discal Cell: A large, distinct, and typically elongate hexagonal discal cell (cell dm) is present in the center of the wing.",
            "Radial Veins: The Radius vein branches into R1, R2+3, R4, and R5. The vein R2+3 is often short and curves forward to meet the Costa.",
            "M-Cu Junction: The base of the cubital fork (where CuA1 and CuA2 diverge) is located near the crossvein m-cu."
        ],
        "Trophic Level": "Predator. Robber flies are voracious predators as adults, capturing other insects in mid-air. Their larvae are also predatory, typically living in soil or rotting wood and feeding on other insect larvae."
    }

    # --- Wing C Analysis ---
    wing_c_description = {
        "Identification": "Family Acrididae (Short-horned Grasshoppers), Order Orthoptera",
        "Venation (Comstock-Needham System)": [
            "Wing Type: This is a tegmen, a leathery and hardened forewing that protects the delicate, membranous hindwing used for flight.",
            "Shape: The wing is long and narrow.",
            "Veins: The main longitudinal veins (C, Sc, R, M, Cu, A) are strong and mostly run parallel to each other.",
            "Crossveins: Numerous crossveins connect the main veins, creating a reticulated (net-like) pattern of small, squarish or irregular cells.",
            "Cells: There are no large, distinct cells like the discal cell in wing B or the specific submarginal cells in wing A."
        ],
        "Trophic Level": "Herbivore. Grasshoppers in the family Acrididae are classic herbivores, feeding on the foliage of various plants."
    }

    # --- Print a summary and conclusion ---
    print("--- Analysis of Insect Wings ---")

    print("\n[Wing A]")
    print(f"Identification: {wing_a_description['Identification']}")
    print("Venation Description:")
    for line in wing_a_description['Venation (Comstock-Needham System)']:
        print(f"  - {line}")
    print(f"Trophic Level: {wing_a_description['Trophic Level']}")

    print("\n[Wing B]")
    print(f"Identification: {wing_b_description['Identification']}")
    print("Venation Description:")
    for line in wing_b_description['Venation (Comstock-Needham System)']:
        print(f"  - {line}")
    print(f"Trophic Level: {wing_b_description['Trophic Level']}")

    print("\n[Wing C]")
    print(f"Identification: {wing_c_description['Identification']}")
    print("Venation Description:")
    for line in wing_c_description['Venation (Comstock-Needham System)']:
        print(f"  - {line}")
    print(f"Trophic Level: {wing_c_description['Trophic Level']}")

    print("\n--- Conclusion ---")
    final_a = wing_a_description['Trophic Level']
    final_b = wing_b_description['Trophic Level']
    final_c = wing_c_description['Trophic Level']

    print(f"A: {final_a}, B: {final_b}, C: {final_c}")
    print("This corresponds to Answer Choice E.")


if __name__ == "__main__":
    analyze_insect_wings()
