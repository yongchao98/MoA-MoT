import re

def check_diels_alder_answer():
    """
    This function verifies the correctness of the provided answer by applying
    the chemical logic and IUPAC name mapping described in the explanation.

    The explanation states the major product is the 'endo-syn' adduct.
    It provides the following mapping to the IUPAC names:
    - 'endo' configuration corresponds to the '(...,4R,7S,...)' descriptor.
    - 'syn' configuration corresponds to the '(...,8s)' descriptor.

    The code will filter the options based on these two rules and check if the
    result matches the given answer.
    """
    # The final answer provided by the LLM
    llm_answer = "C"

    # The options from the question
    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # --- Step 1: Apply the 'endo' rule ---
    # According to the explanation, 'endo' products contain "(4R,7S,".
    endo_candidates = {key for key, name in options.items() if "(4R,7S," in name}
    
    # --- Step 2: Apply the 'syn' rule ---
    # According to the explanation, 'syn' products contain ",8s)".
    syn_candidates = {key for key, name in options.items() if ",8s)" in name}

    # --- Step 3: Find the major product by intersecting the sets ---
    # The major product must be both 'endo' and 'syn'.
    major_product_keys = endo_candidates.intersection(syn_candidates)

    # --- Step 4: Validate the result ---
    # Check if the logic uniquely identifies one product.
    if len(major_product_keys) != 1:
        return (f"Incorrect. The reasoning provided (endo=(4R,7S) and syn=(8s)) "
                f"does not uniquely identify a single product. "
                f"Identified products: {list(major_product_keys)}.")

    # A unique product was found. Check if it matches the LLM's answer.
    predicted_answer = major_product_keys.pop()
    
    if predicted_answer == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The reasoning provided is inconsistent with the final answer. "
                f"The logic points to option {predicted_answer}, but the answer given was {llm_answer}.")

# Execute the check and print the result
result = check_diels_alder_answer()
print(result)