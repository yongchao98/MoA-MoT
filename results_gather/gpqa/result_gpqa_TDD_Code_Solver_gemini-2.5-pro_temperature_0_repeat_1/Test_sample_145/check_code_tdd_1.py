import re

def check_diels_alder_product():
    """
    Checks the correctness of the LLM's answer by applying the known stereochemical
    rules for the Diels-Alder reaction between 5-fluorocyclopenta-1,3-diene and
    maleic anhydride.
    """
    # The options provided in the question
    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # The final answer provided by the LLM
    llm_answer = "A"

    # --- Define Stereochemical Rules ---

    # Rule 1: Endo selectivity is kinetically favored.
    # The endo adduct is identified by the "(...,4S,7R,...)" stereochemical pattern.
    endo_descriptor = "(4S,7R,"

    # Rule 2: Syn selectivity is electronically favored due to the C5-Fluorine.
    # The syn adduct is identified by the "...,8s)" stereochemical pattern.
    syn_descriptor = ",8s)"

    # --- Apply Rules to Find the Correct Product ---
    
    # Find all endo candidates
    endo_candidates = {key for key, name in options.items() if endo_descriptor in name}
    
    # Find all syn candidates
    syn_candidates = {key for key, name in options.items() if syn_descriptor in name}

    # The major product is the intersection of both sets (endo AND syn)
    major_product_keys = endo_candidates.intersection(syn_candidates)

    # --- Validate the Result ---

    # There should be exactly one product that satisfies both rules.
    if len(major_product_keys) != 1:
        return (f"Constraint check failed: The chemical principles (endo and syn selectivity) "
                f"should uniquely identify one major product. However, they identified "
                f"{len(major_product_keys)} products: {list(major_product_keys)}. "
                f"The question options or the principles might be flawed.")

    derived_correct_answer = major_product_keys.pop()

    if derived_correct_answer == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The reasoning based on chemical principles leads to a different answer. "
                f"The major product should be the 'endo-syn' adduct. "
                f"Based on the IUPAC names, the 'endo' options are {endo_candidates} and the 'syn' options are {syn_candidates}. "
                f"The only option that is both is '{derived_correct_answer}'. "
                f"The provided answer was '{llm_answer}'.")

# Execute the check and print the result
result = check_diels_alder_product()
print(result)