import textwrap

def solve_invasive_species_question():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = {
        'A. Apis mellifera (Western honey bee)': {
            'status': 'Introduced and widespread.',
            'impact': 'While essential for agricultural pollination, it is a non-native species that competes with native pollinators for nectar and pollen resources, potentially leading to declines in native bee populations. Its overall impact is complex and debated.'
        },
        'B. Aedes aegypti (Yellow fever mosquito)': {
            'status': 'Introduced and established, particularly in southern New Mexico.',
            'impact': 'This species is a primary vector for serious human diseases, including Zika, dengue, chikungunya, and yellow fever. Its presence poses a significant public health risk, which constitutes a major negative impact on the human population within the ecosystem.'
        },
        'C. Lycorma delicatula (Spotted lanternfly)': {
            'status': 'Not currently established in New Mexico.',
            'impact': 'A highly destructive agricultural and forest pest in the eastern U.S., but it has not yet had an impact in New Mexico as it is not present there.'
        },
        'D. Bombus pascuorum (Common carder bee)': {
            'status': 'Not present in North America.',
            'impact': 'This is a European bumblebee and is not an invasive species in New Mexico.'
        },
        'E. Leptinotarsa decemlineata (Colorado potato beetle)': {
            'status': 'Native to the region (southwestern North America).',
            'impact': 'Although a major agricultural pest of potatoes, it is not an "introduced" invasive species in New Mexico. It spread from its native range as potato agriculture expanded.'
        },
        'F. Maruca vitrata (Bean pod borer)': {
            'status': 'Present as a pest in warmer regions.',
            'impact': 'This moth is an agricultural pest on legumes. While it causes economic damage, its impact in New Mexico is not considered as severe or widespread as the public health threats posed by other species.'
        }
    }

    print("Evaluating the invasive species impact in New Mexico:\n")
    winner = None
    max_impact_reason = ""

    for species, data in species_data.items():
        print(f"Species: {species}")
        print(f" - Status in NM: {data['status']}")
        print(" - Negative Impact:")
        # Use textwrap for clean printing of long impact descriptions
        for line in textwrap.wrap(data['impact'], width=80, initial_indent="   ", subsequent_indent="   "):
            print(line)
        print("-" * 20)

    print("\nConclusion:")
    print("Comparing the valid candidates (introduced and present in New Mexico):")
    print("- Apis mellifera has a negative impact on native pollinators, but this is often offset by its agricultural benefits.")
    print("- Aedes aegypti poses a direct and severe threat to human health by vectoring multiple dangerous diseases.")
    print("\nThe species with the largest negative impact is Aedes aegypti due to its role as a vector for significant human diseases like dengue and Zika in the state.")

solve_invasive_species_question()
<<<B>>>