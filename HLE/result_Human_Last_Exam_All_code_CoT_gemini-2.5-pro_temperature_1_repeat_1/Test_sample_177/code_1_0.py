def solve_immunology_question():
    """
    This function analyzes the roles of various immune receptors to determine
    which one would enable a T cell to act as an antigen-presenting cell (APC).
    """

    # Define the functions of the molecules in the answer choices.
    options = {
        'A': {"name": "CD86", "function": "Co-stimulatory molecule for T cell activation."},
        'B': {"name": "CD80", "function": "Co-stimulatory molecule for T cell activation."},
        'C': {"name": "MHC class I", "function": "Presents endogenous (internal) antigens. Present on almost all nucleated cells, including T cells."},
        'D': {"name": "TIM-4", "function": "Receptor involved in engulfing apoptotic cells, a form of antigen uptake, not presentation."},
        'E': {"name": "MHC class II", "function": "Presents exogenous (external) antigens. The hallmark of professional APCs."}
    }

    # Define the target function we want to engineer into the T cell.
    target_function = "Act as an antigen-presenting cell by presenting exogenous antigens."

    print("Analyzing the requirements for a T cell to act as an Antigen-Presenting Cell (APC):")
    print("---------------------------------------------------------------------------------")
    print(f"1. The primary goal is to enable the T cell to: {target_function}")
    print("2. T cells can already present internal antigens via MHC class I, but to act as a professional APC, they must present external antigens.")

    print("\nEvaluating the options:")
    print("-----------------------")

    correct_option = None
    for key, value in options.items():
        is_match = "exogenous" in value["function"]
        print(f"Option {key} ({value['name']}): {value['function']}")
        if is_match:
            correct_option = key
            print(f"   -> This function matches the key requirement for an APC.")
        else:
            print(f"   -> This is not the primary molecule for presenting external antigens.")

    print("\nConclusion:")
    print("-----------")
    print("To enable a T cell to present antigens it has taken up from the environment (exogenous antigens),")
    print("it must be engineered to express the molecule specifically responsible for that task.")
    print(f"Based on the analysis, {options[correct_option]['name']} is the receptor that presents exogenous antigens, which is the defining characteristic of a professional APC.")
    print(f"Therefore, expressing MHC class II would allow a T cell to function as an APC.")

# Execute the analysis
solve_immunology_question()