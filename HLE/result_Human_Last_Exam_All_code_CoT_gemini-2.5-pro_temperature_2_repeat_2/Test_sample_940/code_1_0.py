import sys

def analyze_invasive_species():
    """
    Analyzes a predefined list of species to determine which has the largest
    negative impact as an invasive in New Mexico.
    """
    species_data = {
        'A': {
            'species_name': 'Apis mellifera',
            'common_name': 'European Honey Bee',
            'status_nm': 'Introduced',
            'impact_level': 2, # Negative impact is debated and less severe than others
            'impact_summary': 'Introduced for agriculture. Can outcompete native pollinators, but is not considered a primary driver of ecosystem collapse. Its negative impacts are generally viewed as less severe than major disease vectors.'
        },
        'B': {
            'species_name': 'Aedes aegypti',
            'common_name': 'Yellow Fever Mosquito',
            'status_nm': 'Invasive',
            'impact_level': 4, # Highest impact due to public health threat
            'impact_summary': 'An established invasive species in southern New Mexico. It is a vector for serious human diseases like Zika, dengue, and chikungunya. Its presence has a major negative public health and economic impact.'
        },
        'C': {
            'species_name': 'Lycorma delicatula',
            'common_name': 'Spotted Lanternfly',
            'status_nm': 'Not Established in NM',
            'impact_level': 0,
            'impact_summary': 'A highly destructive invasive in the eastern US, but it is not currently established in New Mexico. Therefore, its direct impact in the state is nonexistent at this time.'
        },
        'D': {
            'species_name': 'Bombus pascuorum',
            'common_name': 'Common Carder Bee',
            'status_nm': 'Not a Major Invasive in NM',
            'impact_level': 0,
            'impact_summary': 'A European bumblebee that is not considered a significant invasive species in North America, and particularly not in New Mexico.'
        },
        'E': {
            'species_name': 'Leptinotarsa decemlineata',
            'common_name': 'Colorado Potato Beetle',
            'status_nm': 'Native',
            'impact_level': 1, # Negative impact, but not as an introduced species
            'impact_summary': 'Despite being a major agricultural pest, this species is native to a region that includes parts of New Mexico and Colorado. As a native species, it is not classified as an "introduced invasive" here.'
        },
        'F': {
            'species_name': 'Maruca vitrata',
            'common_name': 'Bean Pod Borer',
            'status_nm': 'Minor Pest',
            'impact_level': 1,
            'impact_summary': 'A pest of legume crops in tropical/subtropical areas. While it can cause agricultural damage, it is not cited as having a large-scale negative ecosystem impact in New Mexico compared to other major invasives.'
        }
    }

    print("Analyzing the impact of each species in New Mexico:")
    print("-" * 50)

    highest_impact_level = -1
    worst_offender = None

    for key, data in species_data.items():
        print(f"[{key}] Species: {data['species_name']} ({data['common_name']})")
        print(f"    Status in New Mexico: {data['status_nm']}")
        print(f"    Impact Summary: {data['impact_summary']}\n")

        if data['impact_level'] > highest_impact_level:
            highest_impact_level = data['impact_level']
            worst_offender = key

    print("-" * 50)
    print("Conclusion:")
    if worst_offender:
        winner_data = species_data[worst_offender]
        print(f"The species with the largest negative impact as an introduced invasive in New Mexico is:")
        print(f"({worst_offender}) {winner_data['species_name']}, the {winner_data['common_name']}.")
        print("This is primarily due to its role as a vector for significant human diseases within the state.")
    else:
        print("Could not determine the species with the highest impact.")

if __name__ == '__main__':
    analyze_invasive_species()