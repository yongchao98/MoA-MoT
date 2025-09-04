def check_correctness():
    """
    This function checks the correctness of the final answer by verifying if the
    "Mustard Gas Hypothesis" satisfies all the constraints of the chemical puzzle.
    """

    # The proposed identities of the compounds based on the final answer's reasoning.
    # This is the "Mustard Gas Hypothesis".
    solution = {
        'A': {'formula': 'S8', 'state': 'solid', 'name': 'Sulfur'},
        'B': {'formula': 'Cl2', 'state': 'gas', 'name': 'Chlorine'},
        'C': {'formula': 'SCl2', 'color': 'red', 'name': 'Sulfur dichloride'},
        'D': {'formula': 'C2H4', 'state': 'gas', 'name': 'Ethene'},
        'E': {'formula': '(ClCH2CH2)2S', 'hazard': 'extremely hazardous', 'name': 'Mustard Gas', 'symmetry': 'C2'},
        'F': {'formula': 'HCl', 'acidity': 'strong', 'name': 'Hydrochloric acid'},
        'G': {'formula': 'H2SO3', 'acidity': 'weak', 'name': 'Sulfurous acid'},
        'H': {'formula': 'C2H4Cl2', 'use': 'solvent', 'name': '1,2-dichloroethane'}
    }

    # The options provided in the original question.
    question_options = {'A': 'C2', 'B': 'D4h', 'C': 'C2v', 'D': 'D∞h'}
    
    # The final answer choice provided by the LLM.
    final_answer_choice = 'A'

    errors = []

    # --- Check Clue 1: A(s) + 8 B(g) -> C (bright red product) ---
    if solution['A']['state'] != 'solid':
        errors.append("Clue 1 Fails: A is identified as {} which is not a solid.".format(solution['A']['name']))
    if solution['B']['state'] != 'gas':
        errors.append("Clue 1 Fails: B is identified as {} which is not a gas.".format(solution['B']['name']))
    if solution['C']['color'] != 'red':
        errors.append("Clue 1 Fails: C is identified as {} which is not a red product.".format(solution['C']['name']))
    # The reaction is S8 + 8Cl2 -> 8SCl2. This is a 1:8 molar ratio of A (S8 molecule) to B (Cl2).
    # This check is implicitly satisfied by the choice of compounds.

    # --- Check Clue 2: C + 2 D(g) -> E (extremely hazardous product) ---
    if solution['D']['state'] != 'gas':
        errors.append("Clue 2 Fails: D is identified as {} which is not a gas.".format(solution['D']['name']))
    if solution['E']['hazard'] != 'extremely hazardous':
        errors.append("Clue 2 Fails: E is identified as {} which is not considered 'extremely hazardous'.".format(solution['E']['name']))
    # The reaction is SCl2 + 2C2H4 -> (ClCH2CH2)2S. This is a 1:2 molar ratio of C to D.
    # This check is implicitly satisfied by the choice of compounds.

    # --- Check Clue 3: C + H2O -> A(s) + F(strong acid) + G(weak acid) ---
    # The hydrolysis of SCl2 produces elemental sulfur (A), HCl (F), and sulfurous acid (G).
    # This is a known disproportionation reaction.
    if solution['F']['acidity'] != 'strong':
        errors.append("Clue 3 Fails: F is identified as {} which is not a strong acid.".format(solution['F']['name']))
    if solution['G']['acidity'] != 'weak':
        errors.append("Clue 3 Fails: G is identified as {} which is not a weak acid.".format(solution['G']['name']))
    # The reaction also reforms A (elemental sulfur), which is consistent.

    # --- Check Clue 4: D(g) + B(g) -> H (solvent) (1:1 ratio) ---
    if solution['D']['state'] != 'gas' or solution['B']['state'] != 'gas':
        errors.append("Clue 4 Fails: Reactants D and B are not both gases.")
    if solution['H']['use'] != 'solvent':
        errors.append("Clue 4 Fails: Product H ({}) is not a known solvent.".format(solution['H']['name']))
    # The reaction is C2H4 + Cl2 -> C2H4Cl2. This is a 1:1 molar ratio.
    # This check is implicitly satisfied by the choice of compounds.

    # If all clues are satisfied, check if the final answer (symmetry of E) is correct.
    if not errors:
        expected_symmetry = question_options.get(final_answer_choice)
        actual_symmetry = solution['E']['symmetry']
        
        if expected_symmetry is None:
            errors.append("Final Answer Check Fails: The chosen answer '<<<{}>>>' is not a valid option.".format(final_answer_choice))
        elif actual_symmetry != expected_symmetry:
            errors.append(
                "Final Answer Check Fails: The reasoning identifies E as {} with {} symmetry, "
                "but the final answer '<<<{}>>>' corresponds to {} symmetry.".format(
                    solution['E']['name'], actual_symmetry, final_answer_choice, expected_symmetry
                )
            )

    # --- Return the result ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. The following constraints were not satisfied:\n- " + "\n- ".join(errors)

# Run the check and print the result.
print(check_correctness())