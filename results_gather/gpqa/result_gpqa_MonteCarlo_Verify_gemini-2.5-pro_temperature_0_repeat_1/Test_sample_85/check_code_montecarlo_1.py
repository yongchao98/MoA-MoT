def check_answer(llm_answer_choice):
    """
    Checks the correctness of the selected option for the given chemistry question.

    The function encodes the following chemical principles:
    1.  LiBH4 selectively reduces esters over carboxylic acids. The resulting hydroxy-acid
        cyclizes to a lactone, retaining stereochemistry.
    2.  BH3 selectively reduces carboxylic acids over esters. The resulting hydroxy-ester
        cyclizes to a lactone, retaining stereochemistry.
    """

    # Define the required stereochemistry for the starting materials based on the products
    required_stereochem_A = 'R'  # To get (R)-lactone, A must be (R)
    required_stereochem_B = 'S'  # To get (S)-lactone, B must be (S)

    # Define the stereochemistry proposed by each option for A and B
    options = {
        'A': {'A': 'S', 'B': 'S'},
        'B': {'A': 'R', 'B': 'R'},
        'C': {'A': 'R', 'B': 'S'},
        'D': {'A': 'S', 'B': 'R'}
    }

    if llm_answer_choice not in options:
        return f"Invalid option '{llm_answer_choice}'. The answer must be one of {list(options.keys())}."

    proposed_option = options[llm_answer_choice]
    proposed_stereochem_A = proposed_option['A']
    proposed_stereochem_B = proposed_option['B']

    # Check Reaction A
    # A + LiBH4 -> (R)-lactone
    # This reaction retains stereochemistry, so A must be (R).
    if proposed_stereochem_A != required_stereochem_A:
        return (f"Incorrect for Reaction A. The product is (R)-lactone. "
                f"The reaction with LiBH4 reduces the ester group, and the subsequent lactonization retains the original stereochemistry. "
                f"Therefore, starting material A must have (R) configuration. "
                f"The chosen answer proposes A has ({proposed_stereochem_A}) configuration.")

    # Check Reaction B
    # B + BH3 -> (S)-lactone
    # This reaction retains stereochemistry, so B must be (S).
    if proposed_stereochem_B != required_stereochem_B:
        return (f"Incorrect for Reaction B. The product is (S)-lactone. "
                f"The reaction with BH3 selectively reduces the carboxylic acid, and the subsequent lactonization retains the original stereochemistry. "
                f"Therefore, starting material B must have (S) configuration. "
                f"The chosen answer proposes B has ({proposed_stereochem_B}) configuration.")

    # If both checks pass, the answer is correct.
    return "Correct"

# The provided answer from the other LLM
llm_answer = "C"

# Run the check
result = check_answer(llm_answer)
print(result)