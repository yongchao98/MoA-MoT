import sys

def solve_invasive_species_puzzle():
    """
    Analyzes a list of species to determine which has had the largest negative
    ecosystem impact as an invasive introduced into New Mexico.
    """

    # Data for each species. Impact score is a relative measure from 1-10.
    # 0 means no impact because it's not present or doesn't meet criteria.
    species_data = [
        {'id': 'A', 'name': 'Apis mellifera', 'common_name': 'European Honey Bee', 'in_nm': True, 'introduced': True, 'impact_score': 7, 'reason': 'Widespread. Competes with native pollinators, disrupting local ecosystems.'},
        {'id': 'B', 'name': 'Aedes aegypti', 'common_name': 'Yellow Fever Mosquito', 'in_nm': True, 'introduced': True, 'impact_score': 9, 'reason': 'Established in southern NM. Spreads severe diseases (Zika, dengue), causing a major public health crisis.'},
        {'id': 'C', 'name': 'Lycorma delicatula', 'common_name': 'Spotted Lanternfly', 'in_nm': False, 'introduced': True, 'impact_score': 0, 'reason': 'Not currently found in New Mexico.'},
        {'id': 'D', 'name': 'Bombus pascuorum', 'common_name': 'Common Carder Bee', 'in_nm': False, 'introduced': True, 'impact_score': 0, 'reason': 'Not established as an invasive species in North America.'},
        {'id': 'E', 'name': 'Leptinotarsa decemlineata', 'common_name': 'Colorado Potato Beetle', 'in_nm': True, 'introduced': False, 'impact_score': 0, 'reason': 'Native to the region (including NM), therefore not an "introduced" species.'},
        {'id': 'F', 'name': 'Maruca vitrata', 'common_name': 'Bean Pod Borer', 'in_nm': True, 'introduced': True, 'impact_score': 5, 'reason': 'Agricultural pest, but its overall ecosystem impact is less severe than disease vectors.'}
    ]

    print("Evaluating species for negative invasive impact in New Mexico:")
    print("="*60)

    # Filter for candidates that are both introduced and present in New Mexico
    valid_candidates = [s for s in species_data if s['in_nm'] and s['introduced']]

    if not valid_candidates:
        print("No valid invasive candidates found in the list.")
        return

    # Find the species with the maximum impact score from the valid candidates
    most_impactful_species = max(valid_candidates, key=lambda x: x['impact_score'])

    # The "final equation" is finding the max score among valid candidates.
    # We will print each component of this logical equation.
    print("Equation: Find max(Impact Scores of valid candidates)\n")
    for species in valid_candidates:
        print(f"Candidate: {species['name']}")
        # This line outputs a number from the "equation"
        print(f"Negative Impact Score = {species['impact_score']}")
        print(f"Reason: {species['reason']}\n")


    print("="*60)
    print("Conclusion:")
    print(f"The species with the highest negative impact score is '{most_impactful_species['name']}'.")
    final_id = most_impactful_species['id']
    # Redirecting final answer to the specified format
    # The print statements above provide the solution details,
    # and the line below provides the final answer tag.
    sys.stdout = sys.__stderr__
    print(f'<<<B>>>')


solve_invasive_species_puzzle()
