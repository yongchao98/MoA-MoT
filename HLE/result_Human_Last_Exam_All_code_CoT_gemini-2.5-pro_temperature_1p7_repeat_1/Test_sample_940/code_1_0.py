import sys

def solve_invasive_species_puzzle():
    """
    Analyzes a list of species to determine which has the largest negative
    ecosystem impact as an introduced invasive in New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European Honey Bee)',
            'is_introduced': True,
            'is_in_nm': True,
            'impact_level': 2, # 1=low, 2=moderate, 3=high, 4=very high
            'notes': "Considered invasive as it competes with native pollinators, but also provides critical agricultural benefits. Its net negative impact is debated."
        },
        'B': {
            'name': 'Aedes aegypti (Yellow Fever Mosquito)',
            'is_introduced': True,
            'is_in_nm': True,
            'impact_level': 4,
            'notes': "A vector for serious human diseases including Zika, dengue, and chikungunya. Public health crises are a massive negative ecosystem impact."
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted Lanternfly)',
            'is_introduced': True,
            'is_in_nm': False,
            'impact_level': 4,
            'notes': "A highly destructive pest, but it is not currently established in New Mexico. Therefore, its current impact in the state is zero."
        },
        'D': {
            'name': 'Bombus pascuorum (Common Carder Bee)',
            'is_introduced': True,
            'is_in_nm': False,
            'impact_level': 0,
            'notes': "A European bee that is not established as an invasive species in North America."
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado Potato Beetle)',
            'is_introduced': False,
            'is_in_nm': True,
            'impact_level': 3,
            'notes': "A major agricultural pest, but it is native to the region (SW US/Mexico) and thus not an 'introduced' invasive species in this context."
        },
        'F': {
            'name': 'Maruca vitrata (Bean Pod Borer)',
            'is_introduced': True,
            'is_in_nm': True,
            'impact_level': 2,
            'notes': "An agricultural pest causing economic damage, but its impact is less severe than major public health threats."
        }
    }

    print("Evaluating candidates for the most impactful invasive species in New Mexico:")
    print("="*70)

    valid_candidates = {}
    for key, data in species_data.items():
        print(f"Checking Option {key}: {data['name']}")
        if data['is_introduced'] and data['is_in_nm']:
            print(f"  -> Verdict: QUALIFIED. It is an introduced species present in New Mexico.")
            print(f"  -> Impact Analysis: {data['notes']}")
            valid_candidates[key] = data
        else:
            if not data['is_introduced']:
                print("  -> Verdict: DISQUALIFIED. This species is native to the region.")
            elif not data['is_in_nm']:
                print("  -> Verdict: DISQUALIFIED. This species is not established in New Mexico.")
        print("-" * 70)

    # Determine the winner from the qualified candidates
    if not valid_candidates:
        print("\nNo valid candidates found based on the criteria.")
        return

    winner_key = max(valid_candidates, key=lambda k: valid_candidates[k]['impact_level'])
    winner_data = valid_candidates[winner_key]

    print("\nFinal Conclusion:")
    print(f"Comparing the qualified candidates, '{winner_data['name']}' has the highest impact level.")
    print("Its role as a vector for dangerous human diseases constitutes the largest negative ecosystem impact among the choices.")
    print(f"\nThe correct answer choice is: {winner_key}")


solve_invasive_species_puzzle()