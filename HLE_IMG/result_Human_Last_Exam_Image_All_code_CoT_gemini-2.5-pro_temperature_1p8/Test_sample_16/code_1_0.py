import textwrap

def analyze_wings():
    """
    Analyzes three insect wings, describes their venation, and determines the
    trophic level of the family they belong to.
    """
    wings_data = [
        {
            "id": "A",
            "family": "Ichneumonidae (Hymenoptera)",
            "trophic_level": "Parasitoid",
            "venation_description": (
                "The forewing displays venation characteristic of the family Ichneumonidae. "
                "Key features following the Comstock-Needham system include a prominent dark "
                "pterostigma on the anterior margin (costa). The radial cell (marginal cell), "
                "formed by the Radius (R) vein, is large and fully enclosed. Centrally located is a "
                "small, pentagonal second submarginal cell, famously known as the areolet. This cell "
                "is formed by the convergence of branches of the Radius and Media, and its presence is "
                "a key diagnostic feature. Multiple closed cubital and discal cells are formed by "
                "the Media (M) and Cubitus (Cu) veins and their connecting crossveins."
            )
        },
        {
            "id": "B",
            "family": "Sphecidae or Crabronidae (Hymenoptera)",
            "trophic_level": "Predator",
            "venation_description": (
                "This hymenopteran forewing is characteristic of a hunting wasp. "
                "The venation shows a very long and narrow radial (marginal) cell that runs "
                "along the distal part of the costa. Below the marginal cell, there are three distinct "
                "submarginal cells, formed by branches of the Radius sector (Rs) and connecting crossveins. "
                "The first discoidal cell is also large and well-defined. This venation pattern is typical for "
                "predatory wasps in families like Sphecidae, which require agile flight for hunting."
            )
        },
        {
            "id": "C",
            "family": "Acrididae (Orthoptera)",
            "trophic_level": "Herbivore",
            "venation_description": (
                "This forewing is a modified, leathery structure known as a tegmen, characteristic of the order "
                "Orthoptera (grasshoppers). Its primary role is protective. The venation is simplified compared to a "
                "membranous wing. The main longitudinal veins (Costa, Subcosta, Radius, Media, Cubitus, Anal) are "
                "present but are often thickened and form a straighter, more parallel arrangement. The wing is "
                "coriaceous (leathery) at the base and more membranous with a reticulated pattern towards the apex."
            )
        }
    ]

    print("--- Insect Wing Analysis ---\n")
    for wing in wings_data:
        print(f"Wing ID: {wing['id']}")
        print(f"Identified Family: {wing['family']}")
        print(f"Highest Trophic Level: {wing['trophic_level']}")
        print("Venation Description:")
        # Wrap text for better readability in terminal
        wrapped_text = textwrap.fill(wing['venation_description'], width=80)
        print(f"{wrapped_text}\n")
        print("-" * 30 + "\n")

    # Determine the final answer choice
    final_combination = "A: {wings_data[0]['trophic_level']}, B: {wings_data[1]['trophic_level']}, C: {wings_data[2]['trophic_level']}".format(wings_data=wings_data)
    answer_choices = {
        "A": "A: Herbivore, B: Parasitoid, C: Predator",
        "B": "A: Predator, B: Predator, C: Predator",
        "C": "A: Predator, B: Parasitoid, C: Herbivore",
        "D": "A: Herbivore, B: Predator, C: Parasitoid",
        "E": "A: Parasitoid, B: Predator, C: Herbivore",
        "F": "A: Predator, B: Herbivore, C: Parasitoid",
        "G": "A: Herbivore, B: Predator, C: Herbivore",
        "H": "A: Predator, B: Predator, C: Herbivore",
        "I": "A: Parasitoid, B: Herbivore, C: Predator",
    }
    
    final_answer_letter = ""
    for letter, choice in answer_choices.items():
        if final_combination == choice:
            final_answer_letter = letter
            break

    print("--- Final Conclusion ---")
    print(f"The combined trophic levels are: {final_combination}")
    print(f"This corresponds to answer choice {final_answer_letter}.")


if __name__ == '__main__':
    analyze_wings()
    # The final answer is E
    print("<<<E>>>")