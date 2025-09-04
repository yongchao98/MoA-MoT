def check_organic_chemistry_answer():
    """
    This function checks the correctness of the LLM's answer by codifying the
    chemical principles related to the question.
    """

    # --- Part 1: Analysis of Reactant A ---
    # The reaction is: Cyclohexene + A ---> 8,8-diiodobicyclo[4.2.0]octan-7-one
    # This is a [2+2] cycloaddition between an alkene and a ketene.
    # The product structure indicates the ketene must be diiodoketene (I2C=C=O).
    # The IUPAC name for I2C=C=O is 2,2-diiodoethen-1-one.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # --- Part 2: Analysis of Diene Reactivity Order B ---
    # Diene reactivity in Diels-Alder reactions depends on conformation and electronics.
    # We can assign a reactivity score to each diene based on these principles.
    # Higher score = more reactive.
    dienes_info = {
        # 3: cyclopenta-1,3-diene -> Locked in s-cis, most reactive.
        3: {'score': 4, 'reason': "Locked in the required s-cis conformation."},
        # 1: 2,3-dimethylbuta-1,3-diene -> Two internal EDGs, very reactive.
        1: {'score': 3, 'reason': "Activated by two internal electron-donating groups."},
        # 2: (2E,4E)-hexa-2,4-diene -> Two terminal EDGs, reactive.
        2: {'score': 2, 'reason': "Activated by two terminal electron-donating groups."},
        # 4: (2Z,4Z)-hexa-2,4-diene -> Severe steric hindrance prevents s-cis, unreactive.
        4: {'score': 1, 'reason': "Severe steric hindrance prevents the s-cis conformation."}
    }
    
    # Sort the dienes by their score in descending order to get the correct reactivity sequence.
    correct_reactivity_order_B = sorted(dienes_info.keys(), key=lambda k: dienes_info[k]['score'], reverse=True)

    # --- Part 3: Evaluate the LLM's chosen option ---
    # The LLM chose option B.
    llm_choice = {
        'A': '2,2-diiodoethen-1-one',
        'B': [3, 1, 2, 4]
    }

    # Check if the A part of the answer is correct.
    if llm_choice['A'] != correct_reactant_A:
        return (f"Incorrect. The reactant A is wrong. "
                f"Based on the product of the [2+2] cycloaddition, reactant A must be '{correct_reactant_A}', "
                f"but the answer provided '{llm_choice['A']}'.")

    # Check if the B part of the answer is correct.
    if llm_choice['B'] != correct_reactivity_order_B:
        return (f"Incorrect. The diene reactivity order B is wrong. "
                f"The correct order is {correct_reactivity_order_B} based on conformational and electronic factors. "
                f"The answer provided the order {llm_choice['B']}.")

    # If both parts are correct, the answer is correct.
    return "Correct"

# Run the check
result = check_organic_chemistry_answer()
print(result)