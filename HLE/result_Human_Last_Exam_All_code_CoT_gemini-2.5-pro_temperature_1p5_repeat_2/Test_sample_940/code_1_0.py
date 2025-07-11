import sys

def analyze_invasive_species():
    """
    Analyzes a list of species to determine which has the largest negative
    ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European honey bee)',
            'status_nm': 'Introduced, widespread.',
            'rationale': 'While it competes with native pollinators, it is also crucial for agriculture. Its negative impact is significant but partially offset by benefits.',
            'impact_score': 5
        },
        'B': {
            'name': 'Aedes aegypti (Yellow fever mosquito)',
            'status_nm': 'Introduced, established in southern New Mexico.',
            'rationale': 'A major public health vector for serious diseases like Zika, Dengue, and Chikungunya. This direct threat to human health represents a severe negative impact.',
            'impact_score': 9
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted lanternfly)',
            'status_nm': 'Not present in New Mexico.',
            'rationale': 'This species has not been found in New Mexico, so it has no impact there currently.',
            'impact_score': 0
        },
        'D': {
            'name': 'Bombus pascuorum (Common carder bee)',
            'status_nm': 'Not present in North America.',
            'rationale': 'A European species not found in New Mexico.',
            'impact_score': 0
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado potato beetle)',
            'status_nm': 'Native to the region.',
            'rationale': 'This beetle is native to the Southwestern US (including NM/CO) and therefore does not qualify as an "introduced" invasive species.',
            'impact_score': 0
        },
        'F': {
            'name': 'Maruca vitrata (Cowpea pod borer)',
            'status_nm': 'Introduced, agricultural pest.',
            'rationale': 'This species is an agricultural pest, but its impact is less severe and widespread than the public health threat from Aedes aegypti.',
            'impact_score': 4
        }
    }

    print("Analysis of Species Impact in New Mexico:")
    print("-" * 60)

    # Sort species by impact score for a clear ranking
    sorted_species = sorted(species_data.items(), key=lambda item: item[1]['impact_score'], reverse=True)

    highest_impact_key = sorted_species[0][0]
    
    for key, data in sorted_species:
        print(f"Option {key}: {data['name']}")
        print(f"  - Status: {data['status_nm']}")
        print(f"  - Reason for Score: {data['rationale']}")
        print(f"  - Assigned Negative Impact Score: {data['impact_score']}\n")

    print("\nFinal Conclusion:")
    print("Ranking the species by their negative impact score gives the following equation:")
    
    # Construct and print the "final equation"
    equation_parts = []
    for key, data in sorted_species:
        # Get the first two words of the name for brevity in the equation
        short_name = " ".join(data['name'].split()[:2])
        equation_parts.append(f"{short_name} ({data['impact_score']})")
    
    final_equation_str = " > ".join(equation_parts)
    print(final_equation_str)

    winner_data = species_data[highest_impact_key]
    print(f"\nThe species with the largest negative impact is '{winner_data['name']}' ({highest_impact_key}).")

if __name__ == '__main__':
    analyze_invasive_species()