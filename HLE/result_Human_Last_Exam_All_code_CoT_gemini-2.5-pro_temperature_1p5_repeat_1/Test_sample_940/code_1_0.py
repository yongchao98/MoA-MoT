import sys

def analyze_invasive_species():
    """
    Analyzes a list of potential invasive species in New Mexico to identify the one
    with the largest negative ecosystem impact.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera',
            'common_name': 'European Honey Bee',
            'status_nm': 'Introduced and widespread.',
            'impact_summary': 'Complex impact. While it can outcompete native pollinators, it is also essential for agriculture. Its negative impact is significant but debated and often offset by its economic benefits.'
        },
        'B': {
            'name': 'Aedes aegypti',
            'common_name': 'Yellow Fever Mosquito',
            'status_nm': 'Introduced and established, particularly in southern New Mexico.',
            'impact_summary': 'Major negative impact. This species is a vector for serious human diseases including dengue, chikungunya, Zika, and yellow fever. Its presence poses a direct public health threat, leading to costly control measures and health crises.'
        },
        'C': {
            'name': 'Lycorma delicatula',
            'common_name': 'Spotted Lanternfly',
            'status_nm': 'Not currently known to be established in New Mexico.',
            'impact_summary': 'A highly destructive pest of over 70 plant species in areas where it is established (primarily eastern US). Its potential impact is enormous, but its actual impact in NM is currently negligible.'
        },
        'D': {
            'name': 'Bombus pascuorum',
            'common_name': 'Common Carder Bee',
            'status_nm': 'Not native to North America and not considered a major invasive species in New Mexico.',
            'impact_summary': 'This European bee is not a recognized invasive threat in the region. Its impact is minimal to non-existent.'
        },
        'E': {
            'name': 'Leptinotarsa decemlineata',
            'common_name': 'Colorado Potato Beetle',
            'status_nm': 'Native to North America, including the region of New Mexico.',
            'impact_summary': 'While a major agricultural pest, it is not an "introduced" invasive species in New Mexico. It expanded its range and host plant after the introduction of potatoes.'
        },
        'F': {
            'name': 'Maruca vitrata',
            'common_name': 'Bean Pod Borer',
            'status_nm': 'Present as an agricultural pest.',
            'impact_summary': 'A pest of legume crops. Its impact is primarily agricultural and economic, and it is not typically cited as having the single largest ecosystem-level negative impact in New Mexico compared to major disease vectors.'
        }
    }

    print("Analysis of potential invasive species in New Mexico:")
    print("="*60)

    for choice, data in species_data.items():
        print(f"Choice: {choice}")
        print(f"  Species: {data['name']} ({data['common_name']})")
        print(f"  Status in New Mexico: {data['status_nm']}")
        print(f"  Impact: {data['impact_summary']}")
        print("-"*60)
        # Flush the output to ensure it appears in order.
        sys.stdout.flush()

    print("\nConclusion:")
    print("Of the species listed, Aedes aegypti is the only one that is both introduced to New Mexico and causes a severe, direct negative impact through the transmission of dangerous diseases to the human population. This public health crisis represents the most significant negative impact among the choices.")

analyze_invasive_species()
print("\n<<<B>>>")