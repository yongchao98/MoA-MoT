def check_chemistry_answer():
    """
    Checks the correctness of the selected answer for the given chemistry problem.
    The problem has two parts:
    A) Identifying a reactant in a [2+2] cycloaddition.
    B) Ordering the reactivity of dienes in a Diels-Alder reaction.
    """

    # Define the options from the problem
    options = {
        'A': {
            'reactant_A': "2,2-diiodoethen-1-one",
            'reactivity_B': [3, 1, 2, 4]
        },
        'B': {
            'reactant_A': "4,4-diiodocyclobut-2-en-1-one",
            'reactivity_B': [4, 2, 1, 3]
        },
        'C': {
            'reactant_A': "4,4-diiodocyclobut-2-en-1-one",
            'reactivity_B': [3, 1, 2, 4]
        },
        'D': {
            'reactant_A': "2,2-diiodoethen-1-one",
            'reactivity_B': [4, 2, 1, 3]
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer_key = 'A'
    
    # --- Ground Truth based on Chemical Principles ---

    # Part A: The reaction is a [2+2] cycloaddition. Retrosynthesis of the product
    # 8,8-diiodobicyclo[4.2.0]octan-7-one from cyclohexene reveals the other
    # reactant must be diiodoketene (I2C=C=O), whose IUPAC name is 2,2-diiodoethen-1-one.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # Part B: Diene reactivity in Diels-Alder reactions (most to least reactive):
    # 3. cyclopenta-1,3-diene: Locked s-cis conformation (most reactive).
    # 1. 2,3-dimethylbuta-1,3-diene: Strong internal EDGs, easy s-cis adoption.
    # 2. (2E,4E)-hexa-2,4-diene: Terminal EDGs, less activating than (1).
    # 4. (2Z,4Z)-hexa-2,4-diene: Severe steric hindrance prevents s-cis (least reactive).
    correct_reactivity_B = [3, 1, 2, 4]

    # --- Verification ---
    
    selected_option = options.get(llm_answer_key)
    if not selected_option:
        return f"Invalid answer key '{llm_answer_key}'. Must be one of {list(options.keys())}."

    errors = []

    # Check Part A
    if selected_option['reactant_A'] != correct_reactant_A:
        errors.append(
            f"Constraint on Reactant A is not satisfied. "
            f"The correct reactant is '{correct_reactant_A}', but the answer proposes '{selected_option['reactant_A']}'."
        )

    # Check Part B
    if selected_option['reactivity_B'] != correct_reactivity_B:
        errors.append(
            f"Constraint on Reactivity Order B is not satisfied. "
            f"The correct order is {correct_reactivity_B}, but the answer proposes {selected_option['reactivity_B']}."
        )

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check
result = check_chemistry_answer()
print(result)