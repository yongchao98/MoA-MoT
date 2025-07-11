def find_apc_receptor():
    """
    Analyzes immune cell receptors to determine which one would enable a T cell
    to act as an Antigen-Presenting Cell (APC).
    """
    # Define the biological roles of the potential receptors.
    receptor_info = {
        'A': {'name': 'CD86', 'function': 'Co-stimulatory signal, not antigen presentation.'},
        'B': {'name': 'CD80', 'function': 'Co-stimulatory signal, not antigen presentation.'},
        'C': {'name': 'MHC class I', 'function': 'Presents internal antigens to cytotoxic T cells.'},
        'D': {'name': 'TIM-4', 'function': 'Receptor for clearing apoptotic cells.'},
        'E': {'name': 'MHC class II', 'function': 'Presents external antigens to helper T cells, a key APC function.'}
    }

    # The defining function of an APC we want to engineer into a T cell.
    target_function = "Presents external antigens"

    print("Goal: Engineer a T cell to function as an Antigen-Presenting Cell (APC).")
    print(f"The key function required for this is: '{target_function}'.")
    print("\n--- Evaluating Options ---")

    correct_option = None
    for option, info in receptor_info.items():
        # Check if the receptor's function matches the target APC function.
        is_match = target_function in info['function']
        print(f"Option {option} ({info['name']}): {info['function']} -> Match: {is_match}")
        if is_match:
            correct_option = option
            correct_info = info

    if correct_option:
        print("\n--- Conclusion ---")
        print(f"To make a T cell act as an APC, it needs to present external antigens.")
        print(f"The receptor responsible for this is {correct_info['name']}.")
        print(f"Therefore, the correct choice is Option {correct_option}.")
    else:
        print("\nCould not determine the correct receptor based on the criteria.")

# Run the analysis
find_apc_receptor()