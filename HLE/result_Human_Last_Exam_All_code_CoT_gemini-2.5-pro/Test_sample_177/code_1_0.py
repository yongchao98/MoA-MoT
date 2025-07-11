def find_receptor_for_apc_function():
    """
    This function identifies which receptor would enable a T cell to act as
    an Antigen-Presenting Cell (APC) for helper T cells.
    """
    # Define the primary functions of the candidate receptors
    receptor_functions = {
        "CD86": "Co-stimulation of T cells",
        "CD80": "Co-stimulation of T cells",
        "MHC class I": "Presentation of endogenous antigens to cytotoxic T cells",
        "TIM-4": "Recognition of apoptotic cells",
        "MHC class II": "Presentation of exogenous antigens to helper T cells"
    }

    # The key function required for a cell to act as a professional APC
    # is the ability to present external antigens to helper T cells.
    required_function = "Presentation of exogenous antigens to helper T cells"

    # Identify the correct receptor
    correct_receptor = None
    for receptor, function in receptor_functions.items():
        if function == required_function:
            correct_receptor = receptor
            break

    # Print the reasoning and the result
    print("Goal: To make a T cell act as an Antigen-Presenting Cell (APC).")
    print("A key function of an APC is to present processed external antigens to helper T cells.")
    print(f"This specific function is mediated by the '{required_function}' molecule.")
    print(f"\nEvaluating the options:")
    for receptor, function in receptor_functions.items():
        print(f"- {receptor}: Responsible for '{function}'")

    if correct_receptor:
        print(f"\nConclusion: Engineering a T cell to express {correct_receptor} would enable it to present exogenous antigens to helper T cells, thus allowing it to act as an APC.")
        print(f"Therefore, the correct answer is {correct_receptor}.")
    else:
        print("\nConclusion: Could not identify a suitable receptor from the options.")

find_receptor_for_apc_function()