def find_apc_receptor():
    """
    Analyzes which receptor would enable a T cell to act as an Antigen-Presenting Cell (APC).
    """
    # The primary role of a professional APC is to process external (exogenous) antigens
    # and present them to helper T cells. The receptor for this is MHC class II.
    required_function = "Presenting external antigens to helper T cells"
    print(f"Objective: Find a receptor that enables a T cell to perform the key APC function: '{required_function}'.\n")

    # Dictionary of the options and their primary functions.
    options = {
        "A": {"receptor": "CD86", "function": "A co-stimulatory molecule that provides a 'second signal' for T cell activation, but does not present the antigen."},
        "B": {"receptor": "CD80", "function": "A co-stimulatory molecule, similar to CD86, that helps activate T cells but does not present the antigen."},
        "C": {"receptor": "MHC class I", "function": "Presents internal (endogenous) antigens to cytotoxic T cells. Most cells, including T cells, already have this."},
        "D": {"receptor": "TIM-4", "function": "A receptor involved in recognizing and engulfing apoptotic (dying) cells, not in antigen presentation to T cells."},
        "E": {"receptor": "MHC class II", "function": "Presents external (exogenous) antigens to helper T cells. This is the defining feature of professional APCs."}
    }

    correct_option = None
    print("Evaluating the options:")
    print("-----------------------")
    for key, data in options.items():
        print(f"Option {key}: {data['receptor']}")
        print(f"  - Function: {data['function']}")
        # Check if the function matches the required APC role.
        if "Presents external (exogenous) antigens to helper T cells" in data['function']:
            correct_option = key
            print("  - Verdict: This matches the required function for a T cell to act as an APC.")
        else:
            print("  - Verdict: This does not grant the core APC capability.")
        print()

    print("-----------------------")
    print(f"Conclusion: To engineer a T cell to act as an APC, it must be able to present external antigens.")
    print(f"The receptor that performs this function is {options[correct_option]['receptor']}.")
    print(f"Therefore, the correct answer is option {correct_option}.")

find_apc_receptor()