import re

def check_correctness():
    """
    Checks the correctness of the answer to a chemistry question about epoxide ring-opening.

    The code verifies the following constraints:
    1. Regioselectivity: Attack occurs at the less hindered carbon (C6).
    2. Product Constitution: The resulting structure is a 1,2,4,5-tetramethylcyclohexan-1-ol.
    3. Stereoselectivity: The configuration at the attack site (C6) inverts (S -> R), while others are retained.
    4. IUPAC Naming: The stereochemical descriptors are correctly assigned to the final product's numbering.
    """

    # 1. Define the problem's initial conditions and the proposed answer.
    reactant_stereochem = {'C1': 'R', 'C3': 'R', 'C4': 'R', 'C6': 'S'}
    epoxide_carbon_substitution = {'C1': 'quaternary', 'C6': 'tertiary'}
    
    options = {
        "A": "(1R,2R,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol",
        "B": "(1R,4R,5R)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "C": "(1S,4R,5S)-2,2,4,5-tetramethylcyclohexan-1-ol",
        "D": "(1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol"
    }
    
    # The final answer provided by the LLM to be checked.
    llm_answer_key = "A"
    llm_answer_name = options[llm_answer_key]

    # 2. Apply chemical rules to derive the correct product.

    # Constraint 1: Regioselectivity.
    # Organocuprates attack the less sterically hindered carbon.
    # 'tertiary' is less hindered than 'quaternary'.
    if epoxide_carbon_substitution['C6'] == 'tertiary' and epoxide_carbon_substitution['C1'] == 'quaternary':
        attack_site = 'C6'
    else:
        return "Constraint check failed: The problem definition of steric hindrance is incorrect. C6 should be less hindered than C1."

    # Constraint 2: Product Constitution.
    # Attack at C6 yields a 1,2,4,5-substituted cyclohexanol.
    # Attack at C1 would yield a 2,2,4,5-substituted cyclohexanol.
    if attack_site == 'C6':
        expected_base_name = "1,2,4,5-tetramethylcyclohexan-1-ol"
    else:
        expected_base_name = "2,2,4,5-tetramethylcyclohexan-1-ol"

    if expected_base_name not in llm_answer_name:
        return f"Incorrect. The regioselectivity is wrong. Attack at the less hindered C6 should result in a '{expected_base_name}' structure, but the answer has a different substitution pattern."

    # Constraint 3: Stereoselectivity (Sₙ2 mechanism).
    # Inversion occurs at the attack site (C6), retention elsewhere.
    def invert_config(config):
        return 'S' if config == 'R' else 'R'

    # Map original stereocenters to the new IUPAC numbering of the product.
    # New C1 <- Original C1 (retains config)
    # New C2 <- Original C6 (inverts config)
    # New C4 <- Original C4 (retains config)
    # New C5 <- Original C3 (retains config)
    
    derived_product_stereochem = {
        '1': reactant_stereochem['C1'],
        '2': invert_config(reactant_stereochem['C6']),
        '4': reactant_stereochem['C4'],
        '5': reactant_stereochem['C3']
    }

    # 4. Construct the full name of the derived product and compare.
    derived_stereochem_str = f"(1{derived_product_stereochem['1']},2{derived_product_stereochem['2']},4{derived_product_stereochem['4']},5{derived_product_stereochem['5']})"
    derived_product_full_name = f"{derived_stereochem_str}-{expected_base_name}"

    if derived_product_full_name == llm_answer_name:
        return "Correct"
    else:
        # Find the specific point of error for a more detailed message.
        # Use regex to parse the stereochemistry from the LLM's answer string.
        match = re.match(r'\((\d[RS]),(\d[RS]),(\d[RS]),(\d[RS])\)', llm_answer_name)
        if match:
            llm_stereochem = {
                match.group(1)[0]: match.group(1)[1],
                match.group(2)[0]: match.group(2)[1],
                match.group(3)[0]: match.group(3)[1],
                match.group(4)[0]: match.group(4)[1]
            }
            for center in derived_product_stereochem:
                if derived_product_stereochem[center] != llm_stereochem[center]:
                    return (f"Incorrect. The stereochemistry at C{center} is wrong. "
                            f"Based on the reaction mechanism, it should be {derived_product_stereochem[center]}, "
                            f"but the answer states it is {llm_stereochem[center]}.")
        
        return f"Incorrect. The derived product is {derived_product_full_name}, but the provided answer is {llm_answer_name}."

# Execute the check and print the result.
result = check_correctness()
print(result)