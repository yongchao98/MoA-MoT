import textwrap

def find_receptor_for_apc_function():
    """
    This script determines which receptor would enable a T cell to act as an
    Antigen-Presenting Cell (APC).
    """

    # The defining characteristic of a professional APC is the ability to present
    # processed exogenous antigens to helper T cells.
    key_apc_function = "Presentation of exogenous antigens on MHC class II"
    target_molecule = "MHC class II"

    # The options provided.
    options = {
        'A': 'CD86 (co-stimulatory molecule)',
        'B': 'CD80 (co-stimulatory molecule)',
        'C': 'MHC class I (for endogenous antigens)',
        'D': 'TIM-4 (phagocytosis receptor)',
        'E': 'MHC class II (for exogenous antigens)'
    }

    correct_option = None
    explanation = ""

    # Analyze each option to see if it confers the key APC function.
    for key, description in options.items():
        if target_molecule in description:
            correct_option = key
            explanation = textwrap.dedent(f"""
            Step 1: Identify the core function of a professional Antigen-Presenting Cell (APC).
            The primary role of a professional APC (like a dendritic cell or macrophage) is to take up antigens from the environment (exogenous antigens), process them, and present them to helper T cells. This presentation occurs via a specific receptor called MHC class II.

            Step 2: Understand the normal state of a T cell.
            T cells do not normally express MHC class II. They express MHC class I for presenting internal (endogenous) proteins, but this does not allow them to act as professional APCs for activating other helper T cells.

            Step 3: Evaluate the options to find the molecule that confers APC function.
            - A (CD86) and B (CD80) are co-stimulatory molecules that provide a 'second signal' for T cell activation but do not present the antigen.
            - C (MHC class I) is already present on T cells and presents endogenous, not exogenous, antigens.
            - D (TIM-4) is involved in clearing dead cells, not primary antigen presentation to T cells.
            - E (MHC class II) is the specific molecule used by APCs to present exogenous antigens.

            Step 4: Conclude the solution.
            To enable a T cell to act as an APC, it must be engineered to express the molecule responsible for presenting exogenous antigens. Therefore, introducing MHC class II is the correct approach.

            The correct option is {key}, which corresponds to {options[key]}.
            """)
            break

    print(explanation)
    # The final answer format is requested at the end.
    # The script has identified 'E' as the correct option.
    return f"<<<{correct_option}>>>"

# Execute the function and print the final answer
final_answer = find_receptor_for_apc_function()
print(final_answer)
