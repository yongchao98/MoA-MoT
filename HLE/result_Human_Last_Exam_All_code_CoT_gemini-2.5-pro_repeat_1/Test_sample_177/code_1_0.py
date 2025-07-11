def solve_immunology_question():
    """
    This function analyzes the role of different receptors to determine which one
    would enable a T cell to act as an antigen-presenting cell (APC).
    """

    # Dictionary of choices with their primary functions in immunology.
    receptor_functions = {
        'A': {'name': 'CD86', 'function': 'Co-stimulatory molecule for T cell activation, not antigen presentation.'},
        'B': {'name': 'CD80', 'function': 'Co-stimulatory molecule for T cell activation, not antigen presentation.'},
        'C': {'name': 'MHC class I', 'function': 'Presents endogenous antigens to cytotoxic T cells (CD8+).'},
        'D': {'name': 'TIM-4', 'function': 'Receptor for apoptotic cells, involved in clearance, not presentation.'},
        'E': {'name': 'MHC class II', 'function': 'Presents exogenous antigens to helper T cells (CD4+), the key function of an APC.'}
    }

    # The core requirement for a cell to act as an APC is to present external
    # antigens to helper T cells. We search for the molecule that performs this function.
    correct_choice = None
    for key, value in receptor_functions.items():
        if "presents exogenous antigens" in value['function'].lower():
            correct_choice = key
            break

    # Print the reasoning and the final answer.
    print("Step 1: Understand the goal. The goal is to make a T cell act like an Antigen-Presenting Cell (APC).")
    print("Step 2: Define the key function of an APC. A professional APC presents processed external (exogenous) antigens to helper T cells to initiate an immune response.")
    print("Step 3: Evaluate the options based on this function.")
    for key, value in receptor_functions.items():
        print(f"  - Option {key} ({value['name']}): {value['function']}")
    
    print("\nStep 4: Conclude. The only molecule listed that performs the key APC function of presenting exogenous antigens to helper T cells is MHC class II.")
    
    print(f"\nThe correct answer is {correct_choice}, which corresponds to {receptor_functions[correct_choice]['name']}.")

solve_immunology_question()
print("<<<E>>>")