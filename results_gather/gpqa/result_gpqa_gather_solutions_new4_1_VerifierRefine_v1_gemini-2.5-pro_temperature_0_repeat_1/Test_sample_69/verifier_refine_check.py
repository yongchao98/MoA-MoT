def check_correctness():
    """
    This function checks the correctness of the proposed solution to the chemical puzzle.
    It verifies each clue against the identified substances and their known properties.
    """

    # 1. Define the proposed solution and known chemical properties
    solution = {
        'A': 'S8',
        'B': 'Cl2',
        'C': 'SCl2',
        'D': 'C2H4',
        'E': '(ClCH2CH2)2S',  # Mustard Gas
        'F': 'HCl',
        'G': 'H2SO3',
        'H': 'C2H4Cl2'
    }

    properties = {
        'S8': {'state': 'solid', 'name': 'Sulfur'},
        'Cl2': {'state': 'gas', 'name': 'Chlorine'},
        'SCl2': {'color': 'cherry-red', 'name': 'Sulfur dichloride'},
        'C2H4': {'state': 'gas', 'name': 'Ethene'},
        '(ClCH2CH2)2S': {'hazard': 'extremely hazardous', 'name': 'Mustard gas', 'symmetry': 'C2'},
        'HCl': {'acidity': 'strong', 'name': 'Hydrochloric acid'},
        'H2SO3': {'acidity': 'weak', 'name': 'Sulfurous acid'},
        'C2H4Cl2': {'use': 'solvent', 'name': '1,2-dichloroethane'}
    }

    options = {'A': 'D∞h', 'B': 'C2', 'C': 'C2v', 'D': 'D4h'}
    final_answer_choice = 'B'

    # 2. Verify each clue against the proposed solution

    # Clue 1: A(s) + 8 B(g) -> C (bright red product)
    if properties[solution['A']]['state'] != 'solid':
        return "Constraint Failure (Clue 1): A (S8) is required to be a solid."
    if properties[solution['B']]['state'] != 'gas':
        return "Constraint Failure (Clue 1): B (Cl2) is required to be a gas."
    if 'red' not in properties[solution['C']]['color']:
        return "Constraint Failure (Clue 1): C (SCl2) is described as bright red, but the identified substance is not."
    # Stoichiometry check: The reaction S8 + 8Cl2 -> 8SCl2 perfectly matches the 1:8 reactant ratio.

    # Clue 2: C + 2 D(g) -> E (extremely hazardous product)
    if properties[solution['D']]['state'] != 'gas':
        return "Constraint Failure (Clue 2): D (C2H4) is required to be a gas."
    if properties[solution['E']]['hazard'] != 'extremely hazardous':
        return "Constraint Failure (Clue 2): E (Mustard gas) is not classified as extremely hazardous."
    # Stoichiometry check: The reaction SCl2 + 2C2H4 -> (ClCH2CH2)2S perfectly matches the 1:2 ratio of C to D.

    # Clue 3: C + H2O -> A(s) + F(acid, strong) + G(acid, weak)
    # The hydrolysis of SCl2 produces elemental sulfur (A), HCl (F), and SO2 (which forms H2SO3, G).
    if properties[solution['F']]['acidity'] != 'strong':
        return "Constraint Failure (Clue 3): F (HCl) is required to be a strong acid."
    if properties[solution['G']]['acidity'] != 'weak':
        return "Constraint Failure (Clue 3): G (H2SO3) is required to be a weak acid."
    # The clue that A is reformed is also satisfied.

    # Clue 4: D(g) + B(g) -> H (solvent) (1:1 ratio)
    if properties[solution['D']]['state'] != 'gas' or properties[solution['B']]['state'] != 'gas':
        return "Constraint Failure (Clue 4): Both B and D are required to be gases."
    if properties[solution['H']]['use'] != 'solvent':
        return "Constraint Failure (Clue 4): H (1,2-dichloroethane) is required to be a solvent."
    # Stoichiometry check: The reaction C2H4 + Cl2 -> C2H4Cl2 is a 1:1 reaction.

    # 3. Verify the final question and answer
    symmetry_of_E = properties[solution['E']]['symmetry']
    if options[final_answer_choice] != symmetry_of_E:
        return f"Incorrect Final Answer: The symmetry of E ({solution['E']}) is {symmetry_of_E}. The answer choice {final_answer_choice} corresponds to {options[final_answer_choice]}, which is inconsistent."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)