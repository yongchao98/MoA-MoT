def check_diels_alder_product():
    """
    Checks the correctness of the answer for the Diels-Alder reaction between
    5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The code evaluates the options based on two key principles:
    1. Endo-selectivity: The reaction favors the 'endo' product.
    2. Syn-selectivity: The electronegative fluorine favors a 'syn' approach.
    """

    # Define the options with their key stereochemical descriptors
    options = {
        'A': {'stereochem': ('4R', '7S', '8r'), 'endo_exo': 'exo', 'syn_anti': 'anti'},
        'B': {'stereochem': ('4R', '7S', '8s'), 'endo_exo': 'exo', 'syn_anti': 'syn'},
        'C': {'stereochem': ('4S', '7R', '8s'), 'endo_exo': 'endo', 'syn_anti': 'syn'},
        'D': {'stereochem': ('4S', '7R', '8r'), 'endo_exo': 'endo', 'syn_anti': 'anti'}
    }

    # The provided answer from the other LLM
    given_answer = 'C'

    # --- Step 1: Apply the Endo Rule ---
    # The major product must be the 'endo' adduct.
    endo_candidates = {opt: data for opt, data in options.items() if data['endo_exo'] == 'endo'}

    if given_answer not in endo_candidates:
        return (f"Incorrect. The answer {given_answer} is an '{options[given_answer]['endo_exo']}' product. "
                f"The Diels-Alder reaction is governed by the endo rule, which strongly favors the 'endo' product. "
                f"This eliminates options A and B.")

    # --- Step 2: Apply Syn/Anti Selectivity ---
    # Among the endo products, the major one will be the 'syn' adduct due to the
    # directing effect of the electronegative fluorine atom.
    major_product = None
    for opt, data in endo_candidates.items():
        if data['syn_anti'] == 'syn':
            major_product = opt
            break

    # --- Step 3: Final Verification ---
    if given_answer == major_product:
        return "Correct"
    else:
        return (f"Incorrect. The answer {given_answer} is an 'endo' product, which is a correct first step. "
                f"However, it is the '{options[given_answer]['syn_anti']}' adduct. "
                f"The electronegative fluorine substituent favors the 'syn' approach, not the 'anti' one. "
                f"The major product must be the endo, syn adduct, which is option {major_product}.")

# Execute the check and print the result
result = check_diels_alder_product()
print(result)