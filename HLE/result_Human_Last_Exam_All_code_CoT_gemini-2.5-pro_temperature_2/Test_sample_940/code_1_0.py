import sys

def evaluate_invasive_species():
    """
    Evaluates a list of species to determine which has the largest negative
    impact as an introduced invasive in New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European Honey Bee)',
            'is_introduced': True,
            'is_in_nm': True,
            'reasoning': 'Although it can outcompete native pollinators, it is essential for agriculture. Its impact is complex, not purely negative.',
            'impact_score': 6
        },
        'B': {
            'name': 'Aedes aegypti (Yellow Fever Mosquito)',
            'is_introduced': True,
            'is_in_nm': True,
            'reasoning': 'A vector for serious human diseases like Zika, dengue, and chikungunya. Its presence in southern NM is a major public health threat.',
            'impact_score': 9
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted Lanternfly)',
            'is_introduced': True,
            'is_in_nm': False,
            'reasoning': 'A major pest in the eastern US, but it is not currently established in New Mexico. Therefore, its impact there is zero.',
            'impact_score': 0
        },
        'D': {
            'name': 'Bombus pascuorum (Common Carder Bee)',
            'is_introduced': True,
            'is_in_nm': False,
            'reasoning': 'A European bee that is not known as an invasive species of concern in North America or New Mexico.',
            'impact_score': 0
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado Potato Beetle)',
            'is_introduced': False,
            'is_in_nm': True,
            'reasoning': 'This species is a major agricultural pest, but it is native to the region (SW US/Mexico), not an introduced invasive.',
            'impact_score': 0 # Score is 0 for this exercise as it doesn't meet the "introduced" criteria.
        },
        'F': {
            'name': 'Maruca vitrata (Cowpea Pod Borer)',
            'is_introduced': True,
            'is_in_nm': True,
            'reasoning': 'An agricultural pest of legumes, but not cited as a top-tier invasive threat with widespread impact in New Mexico.',
            'impact_score': 4
        }
    }

    print("Evaluating potential invasive species in New Mexico...\n")

    max_impact_score = -1
    worst_offender_key = None

    for key, data in species_data.items():
        name = data['name']
        score = data['impact_score']
        reason = data['reasoning']
        
        # Only consider species that are introduced and present in NM for the final comparison.
        if data['is_introduced'] and data['is_in_nm']:
            print(f"Candidate: ({key}) {name}\nReasoning: {reason}\nAssigned Negative Impact Score: {score}\n")
            if score > max_impact_score:
                max_impact_score = score
                worst_offender_key = key
        else:
            print(f"Disqualified: ({key}) {name}\nReasoning: {reason}\nAssigned Negative Impact Score: {score}\n")
            
    # Final Result
    if worst_offender_key:
        final_winner = species_data[worst_offender_key]
        winner_name = final_winner['name']
        winner_score = final_winner['impact_score']

        print("--------------------------------------------------")
        print("Final Conclusion:")
        print(f"The species with the largest negative impact is {winner_name}.")
        print("This is represented by the final impact score equation:")
        print(f"{winner_name.split('(')[0].strip()} Impact = {winner_score}")

    else:
        print("Could not determine a winner based on the provided criteria.")


if __name__ == '__main__':
    evaluate_invasive_species()