def check_chemistry_answer():
    """
    Checks the correctness of the selected option for the given chemistry problem.

    The function verifies the answer based on the principle of stereochemistry conservation
    in the described reactions.
    """
    # Define the stereochemistry of the products as given in the question.
    product_A_config = 'R'
    product_B_config = 'S'

    # Based on the chemical principle that the chiral center is not involved in the reaction,
    # its stereochemistry is retained. Therefore, the starting material's configuration
    # must match the product's configuration.
    required_A_config = product_A_config
    required_B_config = product_B_config

    # Define the configurations of the starting materials for each option.
    options = {
        "A": {"A": "R", "B": "R"},
        "B": {"A": "S", "B": "R"},
        "C": {"A": "R", "B": "S"},
        "D": {"A": "S", "B": "S"},
    }

    # The answer provided by the LLM.
    llm_answer = "C"

    # Check if the provided answer is a valid option.
    if llm_answer not in options:
        return f"Invalid option: The provided answer '{llm_answer}' is not one of the possible choices (A, B, C, D)."

    # Get the configurations proposed by the chosen answer.
    chosen_option = options[llm_answer]
    proposed_A_config = chosen_option["A"]
    proposed_B_config = chosen_option["B"]

    # Verify the configuration for starting material A.
    if proposed_A_config != required_A_config:
        return (f"Incorrect. The answer '{llm_answer}' is wrong.\n"
                f"Reason: Reaction A produces an (R) product. Since the stereochemistry is conserved, "
                f"the starting material A must also be (R). "
                f"Option '{llm_answer}' incorrectly proposes that starting material A is ({proposed_A_config}).")

    # Verify the configuration for starting material B.
    if proposed_B_config != required_B_config:
        return (f"Incorrect. The answer '{llm_answer}' is wrong.\n"
                f"Reason: Reaction B produces an (S) product. Since the stereochemistry is conserved, "
                f"the starting material B must also be (S). "
                f"Option '{llm_answer}' incorrectly proposes that starting material B is ({proposed_B_config}).")

    # If both configurations match the requirements, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)