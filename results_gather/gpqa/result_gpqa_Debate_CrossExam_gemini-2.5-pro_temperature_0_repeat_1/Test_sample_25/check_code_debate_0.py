import re

def check_correctness(llm_answer):
    """
    Checks the correctness of the LLM's answer for a two-part chemistry question.

    The function verifies two aspects of the answer:
    1. The identity of reactant A in a [2+2] cycloaddition reaction.
    2. The reactivity order of four dienes in a Diels-Alder type reaction.

    Args:
        llm_answer: A string containing the LLM's full response, which must include the final answer
                    in the format <<<X>>> where X is one of A, B, C, or D.

    Returns:
        A string, "Correct" if the answer is correct according to chemical principles,
        or a string explaining the specific reason why the answer is incorrect.
    """

    # Part 1: Define the ground truth based on established chemical principles.

    # For Part A: Identify reactant A
    # The reaction is: Cyclohexene + A ---> 8,8-diiodobicyclo[4.2.0]octan-7-one
    # This is a [2+2] cycloaddition. The product is a fused 6-membered and 4-membered ring system.
    # Cyclohexene provides the 6-carbon ring. Reactant A must provide a 2-carbon unit
    # that forms the rest of the 4-membered ring.
    # The substituents (a ketone group at C7, two iodine atoms at C8) must come from reactant A.
    # This molecular fragment corresponds to diiodoketene (I2C=C=O).
    # The systematic IUPAC name for diiodoketene is 2,2-diiodoethen-1-one.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # For Part B: Determine the reactivity order of dienes for Diels-Alder reactions.
    # Reactivity is governed by:
    # 1. Conformation: The diene must be in the s-cis conformation. Cyclic dienes locked in s-cis are the most reactive.
    #    Steric hindrance that prevents s-cis formation drastically reduces reactivity.
    # 2. Electronic Effects: Electron-donating groups (EDGs) on the diene increase its reactivity. EDGs on the central
    #    carbons (C2, C3) have a stronger activating effect than on the terminal carbons (C1, C4).
    #
    # Analysis of the given dienes:
    # 3. cyclopenta-1,3-diene: Locked in the ideal s-cis conformation. Most reactive.
    # 1. 2,3-dimethylbuta-1,3-diene: Has two EDGs on the central carbons and can easily adopt the s-cis conformation. Very reactive.
    # 2. (2E,4E)-hexa-2,4-diene: Has two EDGs on the terminal carbons and can adopt the s-cis conformation. Reactive.
    # 4. (2Z,4Z)-hexa-2,4-diene: Severe steric hindrance between the terminal methyl groups prevents it from adopting the s-cis conformation. Least reactive.
    #
    # Therefore, the correct order from most to least reactive is: 3 > 1 > 2 > 4
    correct_reactivity_order_B = [3, 1, 2, 4]

    # Part 2: Parse the LLM's selected option from the provided text.

    # Define the options as given in the question
    options = {
        'A': {"A": "4,4-diiodocyclobut-2-en-1-one", "B": [4, 2, 1, 3]},
        'B': {"A": "4,4-diiodocyclobut-2-en-1-one", "B": [3, 1, 2, 4]},
        'C': {"A": "2,2-diiodoethen-1-one", "B": [3, 1, 2, 4]},
        'D': {"A": "2,2-diiodoethen-1-one", "B": [4, 2, 1, 3]}
    }

    try:
        match = re.search(r'<<<([A-D])>>>', llm_answer)
        if not match:
            return "Failure to parse: The final answer was not found in the required format '<<<X>>>'."
        
        selected_option_key = match.group(1)
        
        if selected_option_key not in options:
            return f"Incorrect: The selected option '{selected_option_key}' is not a valid choice (A, B, C, or D)."
            
        selected_answer = options[selected_option_key]
        
    except Exception as e:
        return f"An error occurred while parsing the LLM's answer: {e}"

    # Part 3: Compare the selected answer with the ground truth and provide feedback.

    # Check Part A: Reactant identification
    if selected_answer["A"] != correct_reactant_A:
        return (
            f"Incorrect: The identification of reactant A is wrong.\n"
            f"The answer chose an option where A is '{selected_answer['A']}'.\n"
            f"The correct reactant A is '{correct_reactant_A}'.\n"
            f"Reasoning: The reaction is a [2+2] cycloaddition. The product, 8,8-diiodobicyclo[4.2.0]octan-7-one, "
            f"is formed from cyclohexene and a ketene. The substituents require the ketene to be "
            f"diiodoketene (I2C=C=O), whose systematic IUPAC name is 2,2-diiodoethen-1-one."
        )

    # Check Part B: Diene reactivity order
    if selected_answer["B"] != correct_reactivity_order_B:
        return (
            f"Incorrect: The reactivity order of dienes (B) is wrong.\n"
            f"The answer chose an option with the order: {selected_answer['B']}.\n"
            f"The correct order from most to least reactive is: {correct_reactivity_order_B}.\n"
            f"Reasoning for the correct order:\n"
            f"1. (3) cyclopenta-1,3-diene: Most reactive as it's rigidly locked in the required s-cis conformation.\n"
            f"2. (1) 2,3-dimethylbuta-1,3-diene: Very reactive due to activating electron-donating groups on its central carbons.\n"
            f"3. (2) (2E,4E)-hexa-2,4-diene: Less reactive than (1) as its activating groups are on the terminal carbons.\n"
            f"4. (4) (2Z,4Z)-hexa-2,4-diene: Least reactive due to severe steric hindrance that prevents it from achieving the s-cis conformation."
        )

    # If all checks pass, the answer is correct.
    return "Correct"