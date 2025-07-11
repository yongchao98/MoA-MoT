import collections

def analyze_wings():
    """
    Analyzes the provided insect wings to determine their family and trophic level,
    then identifies the correct answer choice.
    """
    # Define a data structure to hold the analysis for each wing.
    WingInfo = collections.namedtuple('WingInfo', ['id', 'order', 'family_superfamily', 'venation_description', 'trophic_level'])

    wing_data = [
        WingInfo(
            id='A',
            order='Hymenoptera',
            family_superfamily='Ichneumonoidea (e.g., Ichneumonidae)',
            venation_description=(
                "The wing has a prominent pterostigma on the leading edge (Costa). "
                "The radial vein (R) encloses a large marginal cell. Medial (M) and cubital (Cu) "
                "veins form a series of closed cells, including submarginal and discoidal cells. "
                "This complex venation is characteristic of parasitic wasps."
            ),
            trophic_level='Parasitoid'
        ),
        WingInfo(
            id='B',
            order='Diptera',
            family_superfamily='Asilidae',
            venation_description=(
                "The wing features a large, closed discal cell in the center, from which "
                "medial veins (M) originate. The radial (R) veins are distinct. A key feature is "
                "the long, pointed anal cell at the base of the wing. This venation is "
                "characteristic of robber flies."
            ),
            trophic_level='Predator'
        ),
        WingInfo(
            id='C',
            order='Orthoptera',
            family_superfamily='Caelifera (Grasshoppers)',
            venation_description=(
                "This is an elongate, leathery forewing known as a tegmen. The venation consists of "
                "strong, parallel longitudinal veins (R, M, Cu, etc.) connected by many "
                "crossveins, giving it a reticulate appearance. This structure is typical of grasshoppers."
            ),
            trophic_level='Herbivore'
        )
    ]

    print("Detailed Wing Analysis:")
    print("="*50)
    final_answer_map = {}
    for wing in wing_data:
        print(f"Wing ID: {wing.id}")
        print(f"  Order: {wing.order}")
        print(f"  Probable Family/Superfamily: {wing.family_superfamily}")
        print(f"  Venation (Comstock-Needham System): {wing.venation_description}")
        print(f"  Highest Trophic Level: {wing.trophic_level}")
        print("-" * 50)
        final_answer_map[wing.id] = wing.trophic_level

    # Formulate the final answer based on the analysis
    answer_string = f"A: {final_answer_map['A']}, B: {final_answer_map['B']}, C: {final_answer_map['C']}"
    print(f"\nConclusion: The combined trophic levels are '{answer_string}'.")

    # Match with the given choices
    choices = {
        'A': 'A: Herbivore, B: Parasitoid, C: Predator',
        'B': 'A: Predator, B: Predator, C: Predator',
        'C': 'A: Predator, B: Parasitoid, C: Herbivore',
        'D': 'A: Herbivore, B: Predator, C: Parasitoid',
        'E': 'A: Parasitoid, B: Predator, C: Herbivore',
        'F': 'A: Predator, B: Herbivore, C: Parasitoid',
        'G': 'A: Herbivore, B: Predator, C: Herbivore',
        'H': 'A: Predator, B: Predator, C: Herbivore',
        'I': 'A: Parasitoid, B: Herbivore, C: Predator',
        'J': 'Cannot be determined from the provided information.'
    }
    
    final_choice = 'J' # Default to J if no match is found
    for choice, description in choices.items():
        if description == answer_string:
            final_choice = choice
            break

    print(f"This corresponds to answer choice {final_choice}.")
    print(f"<<<{final_choice}>>>")

# Execute the analysis
analyze_wings()