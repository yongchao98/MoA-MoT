def analyze_invasive_species():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = {
        'A': {
            'scientific_name': 'Apis mellifera',
            'common_name': 'Honey Bee',
            'nm_status': 'Present and widespread.',
            'impact': 'Introduced for pollination and honey production. It can compete with native pollinators but is also economically beneficial. Its negative impact is considered less severe than other major invasive species.'
        },
        'B': {
            'scientific_name': 'Aedes aegypti',
            'common_name': 'Yellow Fever Mosquito',
            'nm_status': 'Present and established, particularly in southern counties like Doña Ana.',
            'impact': 'A major public health vector. It transmits serious human diseases including Zika, dengue, and chikungunya. This direct threat to human health constitutes a significant negative impact on the ecosystem.'
        },
        'C': {
            'scientific_name': 'Lycorma delicatula',
            'common_name': 'Spotted Lanternfly',
            'nm_status': 'Not currently established in New Mexico.',
            'impact': 'A highly destructive pest of many plants, but it is not currently causing an impact in New Mexico as it has not been found there.'
        },
        'D': {
            'scientific_name': 'Bombus pascuorum',
            'common_name': 'Common Carder Bee',
            'nm_status': 'Not considered a significant invasive species in North America or New Mexico.',
            'impact': 'This European species is not a known problem in the state.'
        },
        'E': {
            'scientific_name': 'Leptinotarsa decemlineata',
            'common_name': 'Colorado Potato Beetle',
            'nm_status': 'Native to the broader region of the Southwestern US and Mexico.',
            'impact': 'A major agricultural pest, but it is not an "introduced" species to New Mexico; rather, it is a native pest that adapted to a new food source (potatoes).'
        },
        'F': {
            'scientific_name': 'Maruca vitrata',
            'common_name': 'Spotted Pod Borer',
            'nm_status': 'A known agricultural pest in the Americas.',
            'impact': 'It is a pest of legume crops but is not cited as a primary species causing large-scale negative ecosystem or public health impacts in New Mexico compared to other options.'
        }
    }

    print("Evaluating each species for its negative impact in New Mexico:\n")

    # Print the analysis for each option
    for choice, data in species_data.items():
        print(f"Choice {choice}: {data['common_name']} ({data['scientific_name']})")
        print(f"  Status in NM: {data['nm_status']}")
        print(f"  Impact Summary: {data['impact']}\n")

    # Determine and print the final conclusion
    print("-------------------------")
    print("Conclusion:")
    print("-------------------------")
    print("Based on the analysis, Aedes aegypti is the species from the list with the largest current negative impact in New Mexico.")
    print("Its role as an established vector for serious human diseases represents a direct, ongoing, and severe threat to public health, which is a critical component of the ecosystem.")

if __name__ == '__main__':
    analyze_invasive_species()