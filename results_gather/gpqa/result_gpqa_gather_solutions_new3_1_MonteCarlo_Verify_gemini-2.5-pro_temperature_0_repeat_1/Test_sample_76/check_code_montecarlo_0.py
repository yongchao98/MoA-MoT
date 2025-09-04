import re

def check_reaction_A(product_name_A):
    """
    Checks the product of Reaction A.
    A Wittig rearrangement produces a homoallylic alcohol (HO-C-C-C=C).
    We check if the proposed product fits this structural type.
    """
    # Product from options A and B: 4-methyl-1-phenylpent-3-en-1-ol
    # Structure: HO-C(1)H(Ph)-C(2)H2-C(3)H=C(4)(CH3)-CH3
    # The OH is on C1, the double bond starts at C3. This is a homoallylic alcohol.
    if "4-methyl-1-phenylpent-3-en-1-ol" in product_name_A:
        return True, ""
        
    # Product from options C and D: (Z)-2-methyl-5-phenylpent-2-en-1-ol
    # Structure: HO-C(1)H2-C(2)(CH3)=C(3)H-C(4)H2-C(5)H2-Ph
    # The OH is on C1, the double bond starts at C2. This is an allylic alcohol.
    elif "2-methyl-5-phenylpent-2-en-1-ol" in product_name_A:
        return False, "Reaction A is a Wittig rearrangement, which should produce a homoallylic alcohol. The proposed product, '(Z)-2-methyl-5-phenylpent-2-en-1-ol', is an allylic alcohol."
    
    else:
        return False, f"Unknown product for Reaction A: {product_name_A}"

def check_reaction_B(product_name_B):
    """
    Checks the product of Reaction B.
    The Cope rearrangement is an isomerization, so the degree of saturation must be conserved.
    The starting material is 'hexahydro'. The product must also be 'hexahydro'.
    """
    # The starting material is described as a 'hexahydro' derivative.
    if "hexahydro" in product_name_B:
        return True, ""
    elif "tetrahydro" in product_name_B:
        return False, "Reaction B is a Cope rearrangement (an isomerization). The 'hexahydro' starting material must yield a 'hexahydro' product, not a 'tetrahydro' product as proposed."
    else:
        return False, f"Could not determine the saturation level from the product name for Reaction B: {product_name_B}"

def check_correctness():
    """
    Checks the final answer based on the chemical principles of the reactions.
    """
    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'A'

    # Define the options from the question
    options = {
        'A': {
            'A': "4-methyl-1-phenylpent-3-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[cd]indeno[7,1-gh]azulene"
        },
        'B': {
            'A': "4-methyl-1-phenylpent-3-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'C': {
            'A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'D': {
            'A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[cd]indeno[7,1-gh]azulene"
        }
    }

    chosen_products = options.get(llm_answer_choice)
    if not chosen_products:
        return f"Invalid answer choice '{llm_answer_choice}'. Must be one of {list(options.keys())}."

    product_A_name = chosen_products['A']
    product_B_name = chosen_products['B']

    # Check Reaction A's product type
    is_A_correct, reason_A = check_reaction_A(product_A_name)
    if not is_A_correct:
        return f"Incorrect. The product for Reaction A in option {llm_answer_choice} is inconsistent with the reaction mechanism. Reason: {reason_A}"

    # Check Reaction B's stoichiometry
    is_B_correct, reason_B = check_reaction_B(product_B_name)
    if not is_B_correct:
        return f"Incorrect. The product for Reaction B in option {llm_answer_choice} is inconsistent with the reaction mechanism. Reason: {reason_B}"

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
print(check_correctness())