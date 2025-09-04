def check_correctness_of_chemistry_answer(selected_option: str) -> str:
    """
    Checks the correctness of the selected option for the given chemistry problem.

    The function verifies the answer based on two key chemical principles:
    1. For Reaction A, the Michael addition site on the cyclohexanone determines the
       position of the 'oxo' group in the product's IUPAC name.
    2. For Reaction B, the total number of carbons from the reactants determines the
       parent chain length of the product.
    """

    options = {
        'A': {
            'product_A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'product_B': "3-methyl-4-nitrohexanenitrile"
        },
        'B': {
            'product_A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'product_B': "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'C': {
            'product_A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'product_B': "3-methyl-4-nitrohexanenitrile"
        },
        'D': {
            'product_A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'product_B': "2,3-dimethyl-4-nitrobutanenitrile"
        }
    }

    if selected_option.upper() not in options:
        return f"Invalid option '{selected_option}'. The provided answer must be one of {list(options.keys())}."

    chosen_products = options[selected_option.upper()]
    product_A_name = chosen_products['product_A']
    product_B_name = chosen_products['product_B']

    # Constraint 1: Check Reaction A product structure
    # The reaction is a Michael addition. The enolate of 2-ethyl-2,6-dimethylcyclohexan-1-one
    # can only form at the C6 position, as C2 is a quaternary carbon with no alpha-protons.
    # When naming the product as a substituted propanoate, the ring's point of attachment is C1,
    # making the carbonyl group C2. Thus, the name must contain "2-oxocyclohexyl".
    # Any other position, like "4-oxocyclohexyl", is incorrect.
    is_A_correct = "2-oxocyclohexyl" in product_A_name
    if not is_A_correct:
        return (f"Incorrect: The answer '{selected_option}' is wrong.\n"
                f"Reason: Product A's name is '{product_A_name}', which is structurally incorrect. "
                f"The Michael addition mechanism dictates the formation of a '2-oxocyclohexyl' derivative, not a '4-oxocyclohexyl' one.")

    # Constraint 2: Check Reaction B product structure
    # The reaction is a Michael addition between 1-nitropropane (3 carbons) and
    # (E)-but-2-enenitrile (4 carbons). The longest carbon chain in the product that includes
    # the nitrile group has 6 carbons (CN-C-C(CH3)-C(NO2)-C-C).
    # Therefore, the parent name must be "hexanenitrile". A name like "butanenitrile" is incorrect.
    is_B_correct = "hexanenitrile" in product_B_name
    if not is_B_correct:
        return (f"Incorrect: The answer '{selected_option}' is wrong.\n"
                f"Reason: Product B's name is '{product_B_name}', which has an incorrect carbon backbone. "
                f"The reaction should yield a 6-carbon chain, making the parent name 'hexanenitrile', not 'butanenitrile'.")

    # If both products satisfy the constraints, the answer is correct.
    return "Correct"

# The user provided the final answer as 'C'.
# We will run the check on this specific answer.
final_answer_to_check = "C"
result = check_correctness_of_chemistry_answer(final_answer_to_check)
print(result)