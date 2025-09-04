import re

def get_carbon_count(name: str) -> int:
    """
    Parses a simple IUPAC name to estimate the number of carbon atoms.
    This is a simplified parser sufficient for the given problem.
    """
    count = 0
    # Parent chain contribution
    parent_chains = {
        "undec": 11, "dec": 10, "non": 9, "oct": 8, "hept": 7,
        "hex": 6, "pent": 5, "but": 4, "prop": 3, "eth": 2, "meth": 1
    }
    for chain_name, num in parent_chains.items():
        if chain_name in name:
            count += num
            break  # Assume only one parent chain name part

    # Substituent contribution
    substituents = {
        "methyl": 1, "ethyl": 2, "propyl": 3, "butyl": 4
    }
    # Find all substituent parts, e.g., "di-methyl", "tri-ethyl"
    parts = name.split('-')
    for part in parts:
        for sub_name, num in substituents.items():
            if sub_name in part:
                multiplier = 1
                if part.startswith("di"): multiplier = 2
                if part.startswith("tri"): multiplier = 3
                if part.startswith("tetra"): multiplier = 4
                count += num * multiplier
    return count

def check_correctness_of_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The function verifies the answer based on several chemical principles:
    1.  **Isomerism**: The product of a rearrangement must have the same molecular formula as the reactant.
    2.  **Reaction Mechanism**: The product must be consistent with a Cope rearrangement.
    3.  **Product Type**: A Cope rearrangement of a 1,5-diene yields another 1,5-diene.
    4.  **IUPAC Nomenclature**: The final name must match the derived structure.
    """
    try:
        # --- Step 1: Define the problem and the correctly derived answer ---
        reactant_name = "5-butylnona-2,6-diene"
        # The Cope rearrangement mechanism correctly applied yields:
        expected_product_name = "4-ethyl-3-methyldeca-1,5-diene"

        # --- Step 2: Define the options from the original question ---
        options = {
            "A": "5-ethylundeca-2,6-diene",
            "B": "4-ethyl-3-methyldeca-1,5-diene",
            "C": "5-ethyl-4-methyldeca-2,6-diene",
            "D": "5-ethyl-4-methyldeca-2,6-diene"  # Duplicate as in the prompt
        }

        # --- Step 3: Extract the LLM's chosen option ---
        match = re.search(r'<<<([A-D])>>>', llm_answer_text)
        if not match:
            return "Incorrect. The final answer is not in the required format '<<<X>>>' where X is one of A, B, C, or D."
        
        llm_choice = match.group(1)
        
        if llm_choice not in options:
            return f"Incorrect. The chosen option '{llm_choice}' is not a valid choice from the provided list (A, B, C, D)."

        chosen_product_name = options[llm_choice]

        # --- Step 4: Verify the chosen answer against chemical constraints ---

        # Constraint 1: Isomerism (must have the same number of carbons)
        reactant_carbon_count = get_carbon_count(reactant_name)
        product_carbon_count = get_carbon_count(chosen_product_name)
        if reactant_carbon_count != product_carbon_count:
            return (f"Incorrect. The chosen answer '{chosen_product_name}' is not an isomer of the reactant. "
                    f"The reactant has {reactant_carbon_count} carbons, but the chosen product has {product_carbon_count} carbons.")

        # Constraint 2: Product Type (Cope rearrangement gives a 1,5-diene)
        if "1,5-diene" not in chosen_product_name:
            # Extract the diene type from the chosen name for a more informative message
            diene_type = "unknown diene"
            name_parts = chosen_product_name.split('-')
            if len(name_parts) > 1 and "diene" in name_parts[-1]:
                diene_type = name_parts[-2] + "," + name_parts[-1]

            return (f"Incorrect. The reaction is a Cope rearrangement, which characteristically produces a 1,5-diene. "
                    f"The chosen answer '{chosen_product_name}' is a {diene_type}, which is inconsistent with the mechanism.")

        # Constraint 3: Correct Structure (the name must match the derived product)
        if chosen_product_name != expected_product_name:
            return (f"Incorrect. While the chosen product '{chosen_product_name}' is a 1,5-diene, its substitution pattern does not match "
                    f"the one derived from the Cope rearrangement mechanism. The correct product is '{expected_product_name}'.")

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during checking: {e}"
