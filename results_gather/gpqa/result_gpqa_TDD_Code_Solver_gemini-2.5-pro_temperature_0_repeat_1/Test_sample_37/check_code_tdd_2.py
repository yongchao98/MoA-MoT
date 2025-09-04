import re

def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for a Stork enamine alkylation reaction.
    It simulates the reaction steps based on chemical principles to determine the correct product and reagent sequence,
    and then compares this to the given answer.
    """
    # --- 1. Define Chemical Principles and Expected Outcome ---

    # The reaction starts with an iminium salt of pentan-2-one.
    # The parent ketone is pentan-2-one: CH3(C1)-C(=O)(C2)-CH2(C3)-CH2-CH3
    # It has two alpha-carbons: C1 (methyl, less hindered) and C3 (methylene, more hindered).
    
    # Principle 1: Base and Deprotonation
    # LDA is a strong, bulky base. It will deprotonate the less hindered alpha-carbon (C1)
    # to form the kinetic enamine. This means alkylation will occur at C1.
    alkylation_site = "C1"

    # Principle 2: Alkylation
    # The alkylating agent is ethyl iodide (CH3CH2I), so an ethyl group (CH3CH2-) is added.
    
    # Principle 3: Product Structure
    # Adding an ethyl group to C1 of pentan-2-one results in:
    # (CH3CH2)-CH2-C(=O)-CH2CH2CH3
    # This is a 7-carbon chain (heptane).
    # Numbering from either end places the carbonyl group at position 4.
    correct_product_name = "heptan-4-one"

    # Principle 4: Reagent Sequence
    # The reaction requires three distinct steps in order:
    # (i) Base (LDA) to form the nucleophile.
    # (ii) Alkyl halide (CH3CH2I) for the S_N2 reaction.
    # (iii) Acid (H3O+) for hydrolysis to the final ketone.
    correct_sequence_steps = 3

    # --- 2. Get the LLM's Answer ---
    # The LLM chose option D.
    llm_answer_text = "A = (i) LDA, DME (ii) CH3CH2I (iii) H3O+ B = heptan-4-one"

    # --- 3. Validate the LLM's Answer ---
    try:
        # Parse the answer into reagent and product parts
        if 'B =' not in llm_answer_text:
            return "Incorrect: The answer format is invalid. It must specify both reagents (A) and product (B)."
        
        reagents_part, product_part = llm_answer_text.split('B =')
        answer_product = product_part.strip()

        # Check 3.1: Validate the product
        if correct_product_name.lower() not in answer_product.lower():
            return (f"Incorrect: The product is wrong. The reaction should yield {correct_product_name}, "
                    f"but the answer states the product is '{answer_product}'.")

        # Check 3.2: Validate the reagent sequence
        # Find all step markers like (i), (ii), etc.
        reagent_steps = re.findall(r'\((i+)\)', reagents_part)
        
        if len(reagent_steps) != correct_sequence_steps:
            return (f"Incorrect: The reagent sequence is wrong. The reaction requires {correct_sequence_steps} "
                    f"distinct sequential steps, but the answer provides {len(reagent_steps)}.")

        # Check that the key reagents are in the correct steps
        step_details = re.findall(r'\((i+)\)\s*([^()]+)', reagents_part)
        reagent_map = {step[0]: step[1].upper() for step in step_details}

        if 'LDA' not in reagent_map.get('i', ''):
            return "Incorrect: The base (LDA) must be in step (i) to form the enamine."
        if 'CH3CH2I' not in reagent_map.get('ii', ''):
            return "Incorrect: The alkylating agent (CH3CH2I) must be in step (ii) for the alkylation."
        if 'H3O+' not in reagent_map.get('iii', ''):
            return "Incorrect: The acid (H3O+) for hydrolysis must be in the final step (iii)."

    except Exception as e:
        return f"Incorrect: Failed to parse the answer due to an unexpected format. Error: {e}"

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)