import collections

def check_correctness():
    """
    Checks the correctness of the proposed answer to the chemical puzzle.

    The most consistent pathway identifies the compounds as follows:
    A = S₈ (Sulfur)
    B = Cl₂ (Chlorine)
    C = SCl₂ (Sulfur dichloride)
    D = C₂H₄ (Ethene)
    E = (ClCH₂CH₂)₂S (Mustard gas)
    F = HCl (Hydrochloric acid)
    G = H₂SO₃ (Sulfurous acid)
    H = C₂H₄Cl₂ (1,2-Dichloroethane)

    The final answer is that the symmetry of E (Mustard Gas) is C2, which corresponds to option B.
    """

    # Map of the options provided in the question to their point groups
    QUESTION_OPTIONS = {
        'A': 'Dinfh',
        'B': 'C2',
        'C': 'C2v',
        'D': 'D4h'
    }

    # The proposed solution to be checked
    PROPOSED_SOLUTION = {
        'compounds': {
            'A': 'S8',
            'B': 'Cl2',
            'C': 'SCl2',
            'D': 'C2H4',
            'E': '(ClCH2CH2)2S',
            'F': 'HCl',
            'G': 'H2SO3',
            'H': 'C2H4Cl2'
        },
        'final_choice': 'B'
    }

    # A simple database of chemical properties and reaction facts
    CHEMICAL_DATA = {
        'S8': {'state': 'solid', 'name': 'Sulfur'},
        'Cl2': {'state': 'gas', 'name': 'Chlorine'},
        'SCl2': {'state': 'liquid', 'color': 'red', 'name': 'Sulfur dichloride'},
        'C2H4': {'state': 'gas', 'name': 'Ethene'},
        '(ClCH2CH2)2S': {'hazard': 'extremely hazardous', 'symmetry': 'C2', 'name': 'Mustard Gas'},
        'HCl': {'acidity': 'strong', 'name': 'Hydrochloric acid'},
        'H2SO3': {'acidity': 'weak', 'name': 'Sulfurous acid'},
        'C2H4Cl2': {'type': 'solvent', 'name': '1,2-Dichloroethane'},
    }

    # List to store failure reasons
    errors = []

    # Helper function to get compound formula
    def get_formula(compound_letter):
        return PROPOSED_SOLUTION['compounds'][compound_letter]

    # --- Constraint Checks ---

    # Constraint 1: A is a solid, B is a gas, D is a gas
    if not (CHEMICAL_DATA[get_formula('A')]['state'] == 'solid' and
            CHEMICAL_DATA[get_formula('B')]['state'] == 'gas' and
            CHEMICAL_DATA[get_formula('D')]['state'] == 'gas'):
        errors.append("Constraint Failure: The states of A (solid), B (gas), or D (gas) are incorrect.")

    # Constraint 2: C is a bright red product
    if not CHEMICAL_DATA[get_formula('C')]['color'] == 'red':
        errors.append(f"Constraint Failure: Product C ({CHEMICAL_DATA[get_formula('C')]['name']}) is not red.")

    # Constraint 3: E is an extremely hazardous product
    if not CHEMICAL_DATA[get_formula('E')]['hazard'] == 'extremely hazardous':
        errors.append(f"Constraint Failure: Product E ({CHEMICAL_DATA[get_formula('E')]['name']}) is not described as extremely hazardous.")

    # Constraint 4: F is a strong acid, G is a weak acid
    if not (CHEMICAL_DATA[get_formula('F')]['acidity'] == 'strong' and
            CHEMICAL_DATA[get_formula('G')]['acidity'] == 'weak'):
        errors.append("Constraint Failure: The acid strengths of F (strong) and G (weak) are incorrect.")

    # Constraint 5: H is a solvent
    if not CHEMICAL_DATA[get_formula('H')]['type'] == 'solvent':
        errors.append(f"Constraint Failure: Product H ({CHEMICAL_DATA[get_formula('H')]['name']}) is not a solvent.")

    # Constraint 6: Reaction 1 (A + 8B -> C)
    # The reaction is S₈(s) + 8Cl₂(g) -> 8SCl₂(l). This fits the 1:8 stoichiometry perfectly.
    # We will represent this as a check for this specific known reaction.
    if not (get_formula('A') == 'S8' and get_formula('B') == 'Cl2' and get_formula('C') == 'SCl2'):
        errors.append("Constraint Failure: Reaction 1 (A + 8B -> C) is not consistent with S₈ + 8Cl₂ -> 8SCl₂.")

    # Constraint 7: Reaction 2 (C + 2D -> E)
    # The reaction is SCl₂(l) + 2C₂H₄(g) -> (ClCH₂CH₂)₂S(l). This is the Levinstein process.
    if not (get_formula('C') == 'SCl2' and get_formula('D') == 'C2H4' and get_formula('E') == '(ClCH2CH2)2S'):
        errors.append("Constraint Failure: Reaction 2 (C + 2D -> E) is not consistent with SCl₂ + 2C₂H₄ -> Mustard Gas.")

    # Constraint 8: Reaction 3 (C + H2O -> A + F + G)
    # The hydrolysis of SCl₂ produces S, HCl, and SO₂ (which forms H₂SO₃).
    if not (get_formula('C') == 'SCl2' and get_formula('A') == 'S8' and get_formula('F') == 'HCl' and get_formula('G') == 'H2SO3'):
        errors.append("Constraint Failure: Reaction 3 (hydrolysis of C) does not produce A, F, and G as proposed.")

    # Constraint 9: Reaction 4 (D + B -> H, 1:1)
    # The reaction is C₂H₄(g) + Cl₂(g) -> C₂H₄Cl₂(l). This is a 1:1 reaction.
    if not (get_formula('D') == 'C2H4' and get_formula('B') == 'Cl2' and get_formula('H') == 'C2H4Cl2'):
        errors.append("Constraint Failure: Reaction 4 (D + B -> H) is not consistent with C₂H₄ + Cl₂ -> C₂H₄Cl₂.")

    # Final Check: Molecular Symmetry of E
    formula_E = get_formula('E')
    known_symmetry = CHEMICAL_DATA[formula_E]['symmetry']
    proposed_choice = PROPOSED_SOLUTION['final_choice']
    symmetry_from_choice = QUESTION_OPTIONS[proposed_choice]

    # Case-insensitive comparison for point groups
    if known_symmetry.lower() != symmetry_from_choice.lower():
        errors.append(f"Final Answer Failure: The symmetry of E ({formula_E}) is {known_symmetry}, "
                      f"but the chosen answer '{proposed_choice}' corresponds to {symmetry_from_choice}.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "\n".join(["Incorrect. Reason(s):"] + errors)

# Run the check and print the result
result = check_correctness()
print(result)