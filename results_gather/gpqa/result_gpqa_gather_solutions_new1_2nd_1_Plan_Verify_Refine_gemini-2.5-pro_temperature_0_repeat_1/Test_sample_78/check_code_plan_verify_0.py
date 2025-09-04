def check_correctness():
    """
    This function checks the correctness of the provided answer by programmatically
    applying the constraints from the chemistry problem.

    The logic follows these steps:
    1.  Define the properties of the four candidate compounds (A, B, C, D), including
        their molecular formula and whether they contain a p-tolyl group.
    2.  Establish the constraints for the starting material (Compound X) based on the
        problem description:
        -   Constraint 1: The molecular formula must be C11H12O.
        -   Constraint 2: Based on the product's NMR data showing a p-tolyl group,
            the starting material must also contain a p-tolyl group, as the reaction
            is a rearrangement.
    3.  Filter the candidate compounds against these two constraints.
    4.  Compare the single valid option found with the provided final answer.
    """

    # The final answer provided by the LLM analysis to be checked.
    llm_final_answer = 'C'

    # --- Data and Constraints from the Question ---

    # Constraint 1: Required molecular formula for Compound X
    required_formula = {'C': 11, 'H': 12, 'O': 1}

    # Constraint 2: The product has a p-tolyl group, so the starting material must too.
    requires_ptolyl_group = True

    # Properties of the candidate compounds
    options_data = {
        'A': {
            'name': "2-methyl-3-styryloxirane",
            'formula': {'C': 11, 'H': 12, 'O': 1},
            'has_ptolyl_group': False  # "styryl" implies a simple phenyl group
        },
        'B': {
            'name': "2-(1-phenylprop-1-en-2-yl)oxirane",
            'formula': {'C': 11, 'H': 12, 'O': 1},
            'has_ptolyl_group': False  # "phenyl" is explicit
        },
        'C': {
            'name': "2-(4-methylstyryl)oxirane",
            'formula': {'C': 11, 'H': 12, 'O': 1},
            'has_ptolyl_group': True  # "4-methylstyryl" is a p-tolyl derivative
        },
        'D': {
            'name': "2-styrylepoxide",
            'formula': {'C': 10, 'H': 10, 'O': 1}, # Incorrect formula
            'has_ptolyl_group': False  # "styryl" implies a simple phenyl group
        }
    }

    # --- Verification Logic ---

    valid_options = []
    elimination_reasons = {}

    for option_key, data in options_data.items():
        is_valid = True
        reasons = []

        # Check Constraint 1: Molecular Formula
        if data['formula'] != required_formula:
            is_valid = False
            reasons.append(f"has the wrong molecular formula ({data['formula']})")

        # Check Constraint 2: Structural Feature (p-tolyl group)
        if requires_ptolyl_group and not data['has_ptolyl_group']:
            is_valid = False
            reasons.append("lacks the required p-tolyl group")

        if is_valid:
            valid_options.append(option_key)
        else:
            elimination_reasons[option_key] = " and ".join(reasons)

    # --- Final Evaluation ---

    # There should be exactly one valid option
    if len(valid_options) != 1:
        return f"Analysis failed: Found {len(valid_options)} valid options ({', '.join(valid_options)}). Expected exactly one."

    correct_option = valid_options[0]

    if llm_final_answer == correct_option:
        return "Correct"
    else:
        reason = f"Incorrect. The provided answer is {llm_final_answer}, but the correct answer is {correct_option}. "
        reason += f"Option {llm_final_answer} is incorrect because it {elimination_reasons.get(llm_final_answer, 'fails the problem constraints')}. "
        reason += f"Option {correct_option} is the only one that has the correct molecular formula and contains the p-tolyl group necessary for the rearrangement."
        return reason

# Run the check and print the result
result = check_correctness()
print(result)