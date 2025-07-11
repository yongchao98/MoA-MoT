import sys

def analyze_invasive_species_impact():
    """
    Analyzes the negative impact of several species as introduced invasives in New Mexico.

    The function uses a predefined dictionary containing information about each species.
    It assigns a qualitative 'impact_score' to represent the severity of the negative
    impact in New Mexico. A score of 0 is used for species that are not present, native,
    or not considered significant invasives in the state.
    """
    species_data = {
        'A': {
            'scientific_name': 'Apis mellifera',
            'common_name': 'Western honey bee',
            'status_in_nm': 'Introduced and managed. Can compete with native pollinators but is also a key agricultural pollinator.',
            'impact_score': 3
        },
        'B': {
            'scientific_name': 'Aedes aegypti',
            'common_name': 'Yellow fever mosquito',
            'status_in_nm': 'Established invasive in southern NM. A major public health vector for diseases like Zika, dengue, and chikungunya.',
            'impact_score': 9
        },
        'C': {
            'scientific_name': 'Lycorma delicatula',
            'common_name': 'Spotted lanternfly',
            'status_in_nm': 'Not currently established in New Mexico.',
            'impact_score': 0
        },
        'D': {
            'scientific_name': 'Bombus pascuorum',
            'common_name': 'Common carder bee',
            'status_in_nm': 'Not established in North America.',
            'impact_score': 0
        },
        'E': {
            'scientific_name': 'Leptinotarsa decemlineata',
            'common_name': 'Colorado potato beetle',
            'status_in_nm': 'Native to the region (including New Mexico), therefore not an "introduced" invasive species.',
            'impact_score': 0 # Score is 0 in the context of an *introduced* invasive
        },
        'F': {
            'scientific_name': 'Maruca vitrata',
            'common_name': 'Spotted pod borer',
            'status_in_nm': 'Not considered a major invasive pest in New Mexico.',
            'impact_score': 1
        }
    }

    max_impact_score = -1
    most_impactful_species_key = None

    print("Evaluating Negative Impact of Introduced Species in New Mexico")
    print("============================================================")
    print("The 'equation' below compares the impact score of each species.")

    # Iterate through the data to find the max score, printing each 'number' (score)
    for key, data in species_data.items():
        score = data['impact_score']
        print(f"\n- Species ({key}): {data['scientific_name']} ({data['common_name']})")
        print(f"  Status: {data['status_in_nm']}")
        print(f"  Negative Impact Score = {score}") # This is the number for the 'equation'
        if score > max_impact_score:
            max_impact_score = score
            most_impactful_species_key = key

    print("\n============================================================")
    print("Conclusion:")
    if most_impactful_species_key:
        winner = species_data[most_impactful_species_key]
        print(f"The species with the largest negative impact is '{winner['scientific_name']}' with a score of {max_impact_score}.")
        print(f"This corresponds to answer choice {most_impactful_species_key}.")
    else:
        print("Could not determine the most impactful species.")

if __name__ == '__main__':
    analyze_invasive_species_impact()