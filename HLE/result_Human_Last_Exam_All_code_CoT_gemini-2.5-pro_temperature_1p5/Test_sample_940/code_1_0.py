import textwrap

def analyze_invasive_species_impact():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European Honey Bee)',
            'status_nm': 'Established Invasive',
            'impact': ('Widespread and established. Competes directly with New Mexico\'s '
                       'diverse native pollinator species for resources. This competition '
                       'can lead to declines in native bee populations and disrupt natural '
                       'pollination cycles for native plants.')
        },
        'B': {
            'name': 'Aedes aegypti (Yellow Fever Mosquito)',
            'status_nm': 'Established Invasive',
            'impact': ('Established in southern New Mexico. It is a primary vector for '
                       'dangerous human diseases like Zika, dengue, and chikungunya. Its '
                       'impact is a significant public health threat, necessitating costly '
                       'monitoring and control measures. A widespread disease outbreak would '
                       'have a devastating impact on the human population, a key part of '
                       'the ecosystem.')
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted Lanternfly)',
            'status_nm': 'Not Established',
            'impact': ('A major invasive threat in the eastern US, but it is not currently '
                       'established in New Mexico. Therefore, it has not had an impact there.')
        },
        'D': {
            'name': 'Bombus pascuorum (Common Carder Bee)',
            'status_nm': 'Not Present in North America',
            'impact': ('A European species that is not known to be present or invasive in '
                       'North America, including New Mexico.')
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado Potato Beetle)',
            'status_nm': 'Native Species',
            'impact': ('This beetle is native to the southwestern US / Mexico region. '
                       'It became a pest after adapting to potatoes. As it is native, '
                       'it is not an "introduced" invasive species in New Mexico.')
        },
        'F': {
            'name': 'Maruca vitrata (Bean Pod Borer)',
            'status_nm': 'Minor Pest',
            'impact': ('This moth is a pest of legume crops but is not considered a '
                       'primary invasive species causing widespread ecosystem-level '
                       'damage in New Mexico compared to other major invasives.')
        }
    }

    print("Step 1: Evaluating each species in the context of New Mexico.\n")

    valid_candidates = {}
    for key, data in species_data.items():
        if data['status_nm'] in ['Established Invasive']:
            valid_candidates[key] = data

    print("Step 2: Filtering out species that are not established invasives in NM.")
    print("The remaining candidates for largest impact are:\n")
    for key, data in valid_candidates.items():
        print(f"  - {key}: {data['name']}")
    print("\n" + "="*50 + "\n")

    print("Step 3: Comparing the impacts of the valid candidates.\n")

    print(f"Analysis of Candidate A: {valid_candidates['A']['name']}")
    print(textwrap.fill(valid_candidates['A']['impact'], width=80))
    print("\n" + "-"*50 + "\n")
    print(f"Analysis of Candidate B: {valid_candidates['B']['name']}")
    print(textwrap.fill(valid_candidates['B']['impact'], width=80))
    print("\n" + "="*50 + "\n")

    print("Step 4: Conclusion.\n")
    conclusion = ("While the European Honey Bee (A) causes significant harm to native ecosystems by "
                  "outcompeting native pollinators, the Yellow Fever Mosquito (B) poses a direct, "
                  "severe, and potentially lethal threat to human health through disease transmission. "
                  "This significant risk to the human population is considered the largest negative "
                  "impact among the given choices.")
    print(textwrap.fill(conclusion, width=80))
    
    final_answer_key = 'B'
    print(f"\nFinal Conclusion: The species with the largest negative impact is {species_data[final_answer_key]['name']}.")

if __name__ == '__main__':
    analyze_invasive_species_impact()