def check_correctness():
    """
    This function verifies the proposed solution to the chemical puzzle by checking
    each constraint programmatically.
    """

    # The proposed solution to be verified.
    # This is the "Mustard Gas System" hypothesis.
    solution = {
        'A': 'S8',
        'B': 'Cl2',
        'C': 'SCl2',
        'D': 'C2H4',
        'E': '(ClCH2CH2)2S',
        'F': 'HCl',
        'G': 'H2SO3',
        'H': 'C2H4Cl2',
    }

    # A simple database of chemical properties relevant to the puzzle.
    properties = {
        'S8': {'name': 'Sulfur', 'state': 'solid'},
        'Cl2': {'name': 'Chlorine', 'state': 'gas'},
        'SCl2': {'name': 'Sulfur dichloride', 'color': 'cherry-red', 'state': 'liquid'},
        'C2H4': {'name': 'Ethene', 'state': 'gas'},
        '(ClCH2CH2)2S': {'name': 'Mustard gas', 'hazard': 'extremely hazardous', 'point_group': 'C2'},
        'HCl': {'name': 'Hydrochloric acid', 'acidity': 'strong'},
        'H2SO3': {'name': 'Sulfurous acid', 'acidity': 'weak'},
        'C2H4Cl2': {'name': '1,2-dichloroethane', 'use': 'solvent'},
    }
    
    # The final answer from the proposed solution.
    final_answer_choice = 'D'
    final_answer_symmetry = 'C2'

    failures = []

    # --- Check Constraint 1: A(s) + 8 B(g) -> C (bright red product) ---
    # Proposed: S8(s) + 8 Cl2(g) -> 8 SCl2
    if properties[solution['A']]['state'] != 'solid':
        failures.append("Constraint 1: A ({}) is not a solid.".format(solution['A']))
    if properties[solution['B']]['state'] != 'gas':
        failures.append("Constraint 1: B ({}) is not a gas.".format(solution['B']))
    # The stoichiometry 1 mole of S8 to 8 moles of Cl2 matches "8 equivalents".
    if 'red' not in properties[solution['C']]['color']:
        failures.append("Constraint 1: C ({}) is not a 'bright red product'.".format(solution['C']))

    # --- Check Constraint 2: C + 2 D(g) -> E (extremely hazardous product) ---
    # Proposed: SCl2 + 2 C2H4(g) -> (ClCH2CH2)2S
    if properties[solution['D']]['state'] != 'gas':
        failures.append("Constraint 2: D ({}) is not a gas.".format(solution['D']))
    # The stoichiometry 1 mole of SCl2 to 2 moles of C2H4 matches "2 equivalents".
    if properties[solution['E']]['hazard'] != 'extremely hazardous':
        failures.append("Constraint 2: E ({}) is not an 'extremely hazardous product'.".format(solution['E']))

    # --- Check Constraint 3: C + H₂O -> A(s) + F(strong acid) + G(weak acid) ---
    # Proposed: Hydrolysis of SCl2 produces S8(s), HCl(strong), and H2SO3(weak).
    # The reaction 2SCl₂ + 2H₂O → SO₂ + 4HCl + S confirms these products.
    # SO₂ in water gives H₂SO₃.
    if properties[solution['A']]['state'] != 'solid':
        failures.append("Constraint 3: The reformed product A ({}) is not a solid.".format(solution['A']))
    if properties[solution['F']]['acidity'] != 'strong':
        failures.append("Constraint 3: Product F ({}) is not a strong acid.".format(solution['F']))
    if properties[solution['G']]['acidity'] != 'weak':
        failures.append("Constraint 3: Product G ({}) is not a weak acid.".format(solution['G']))

    # --- Check Constraint 4: D(g) + B(g) -> H (solvent) (1:1 ratio) ---
    # Proposed: C2H4(g) + Cl2(g) -> C2H4Cl2
    if properties[solution['D']]['state'] != 'gas':
        failures.append("Constraint 4: D ({}) is not a gas.".format(solution['D']))
    if properties[solution['B']]['state'] != 'gas':
        failures.append("Constraint 4: B ({}) is not a gas.".format(solution['B']))
    # The reaction is 1:1.
    if properties[solution['H']]['use'] != 'solvent':
        failures.append("Constraint 4: Product H ({}) is not a solvent.".format(solution['H']))

    # --- Check the final question: What is the molecular symmetry group of E? ---
    # Proposed: E is Mustard Gas, point group is C2. The answer choice is D.
    point_group = properties[solution['E']]['point_group']
    if point_group != final_answer_symmetry:
        failures.append("Final Question Check: The point group of E ({}) is {}, but the answer claims it is {}.".format(solution['E'], point_group, final_answer_symmetry))
    
    options = {'A': 'D∞h', 'B': 'C2v', 'C': 'D4h', 'D': 'C2'}
    if options.get(final_answer_choice) != final_answer_symmetry:
        failures.append("Final Question Check: The final answer choice '{}' does not correspond to the required symmetry '{}'.".format(final_answer_choice, final_answer_symmetry))

    # --- Final Verdict ---
    if not failures:
        return "Correct"
    else:
        # Return the first failure found
        return "Incorrect: " + failures[0]

# Execute the check
result = check_correctness()
print(result)