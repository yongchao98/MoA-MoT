import re

def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for a chemical reaction question.
    It simulates the chemical logic based on the rules of organic chemistry to derive the correct product
    and compares it to the given answer.
    """

    # --- Problem and Answer Definition ---
    # The question describes the reaction of an epoxide with an organocuprate.
    # Reactant: (1R,3R,4R,6S)-1,3,4-trimethyl-7-oxabicyclo[4.1.0]heptane
    # Reagent: Me2CuLi (provides a methyl nucleophile)
    # The final answer provided by the LLM to be checked is C.
    final_answer_letter = "C"
    options = {
        "A": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "B": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "C": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "D": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol"
    }
    llm_answer_text = options.get(final_answer_letter)

    if not llm_answer_text:
        return f"Invalid answer letter '{final_answer_letter}'. Please choose from A, B, C, or D."

    # --- Chemical Logic Simulation ---

    # Step 1: Analyze Regioselectivity (Site of Attack)
    # The rule is that the nucleophile attacks the less sterically hindered carbon of the epoxide.
    # The epoxide is across C1 and C6.
    # C1 is a quaternary carbon (bonded to a methyl group), which is more hindered.
    # C6 is a tertiary carbon (bonded to a hydrogen), which is less hindered.
    # Therefore, the attack occurs at C6.
    attack_site = 6

    # Step 2: Determine Product Constitution
    # Attack at C6 opens the ring to form an alcohol at C1 and adds a new methyl group at C6.
    # For IUPAC naming, the carbon with the -OH group is the new C1. The original C6 becomes the new C2.
    # The product is a 1,2,4,5-tetramethylcyclohexan-1-ol.
    expected_constitution = "1,2,4,5-tetramethylcyclohexan-1-ol"

    # Constraint Check 1: Constitution
    if expected_constitution not in llm_answer_text:
        return (f"Incorrect. The product's constitution is wrong. "
                f"The regioselectivity rule (attack at less hindered carbon C6) dictates that "
                f"the product should be a '{expected_constitution}', but the answer is '{llm_answer_text}'.")

    # Step 3: Analyze Stereoselectivity
    # The reaction is an Sₙ2 attack, which causes inversion of configuration at the attacked carbon (C6)
    # and retention of configuration at the other chiral centers (C1, C3, C4).
    initial_stereochem = {'C1': 'R', 'C3': 'R', 'C4': 'R', 'C6': 'S'}
    
    # Determine the expected stereochemistry of the product.
    # Note the renumbering: old C1->new C1, old C6->new C2, old C3->new C5, old C4->new C4.
    expected_product_stereochem = {}
    # New C1 (from old C1) retains its configuration.
    expected_product_stereochem['C1'] = initial_stereochem['C1']  # R
    # New C2 (from old C6) inverts its configuration.
    expected_product_stereochem['C2'] = 'R' if initial_stereochem['C6'] == 'S' else 'S'  # S -> R
    # New C4 (from old C4) retains its configuration.
    expected_product_stereochem['C4'] = initial_stereochem['C4']  # R
    # New C5 (from old C3) retains its configuration.
    expected_product_stereochem['C5'] = initial_stereochem['C3']  # R

    # Step 4: Assemble the full expected product name and compare.
    expected_stereochem_str = (f"(1{expected_product_stereochem['C1']},"
                               f"2{expected_product_stereochem['C2']},"
                               f"4{expected_product_stereochem['C4']},"
                               f"5{expected_product_stereochem['C5']})")
    expected_full_name = f"{expected_stereochem_str}-{expected_constitution}"

    # Constraint Check 2: Stereochemistry
    if llm_answer_text == expected_full_name:
        return "Correct"
    else:
        # Parse the stereochemistry from the LLM's answer to provide a detailed error message.
        try:
            answer_stereochem_str = re.search(r'\((.*?)\)', llm_answer_text).group(1)
            answer_configs = {f'C{int(c[0])}': c[1] for c in re.findall(r'(\d+)([RS])', answer_stereochem_str)}
            
            errors = []
            for center, config in expected_product_stereochem.items():
                if answer_configs.get(center) != config:
                    rule = "inverted" if center == 'C2' else "retained"
                    errors.append(f"the configuration at {center} should be {config} ({rule}) but was {answer_configs.get(center, 'not specified')}")
            
            if errors:
                return f"Incorrect. The stereochemistry is wrong. The following constraints were not satisfied: {'; '.join(errors)}."
            else:
                # This case should not be reached if names don't match but stereochem does.
                return f"Incorrect. Expected '{expected_full_name}' but got '{llm_answer_text}'."

        except (AttributeError, IndexError):
            # Fallback error message if parsing fails.
            return f"Incorrect. The expected product is '{expected_full_name}', but the provided answer is '{llm_answer_text}'."

# Execute the check
result = check_chemistry_answer()
print(result)