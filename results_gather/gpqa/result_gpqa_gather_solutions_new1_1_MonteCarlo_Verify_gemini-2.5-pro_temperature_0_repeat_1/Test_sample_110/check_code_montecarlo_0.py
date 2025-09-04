import sys
import io

def check_correctness():
    """
    This function checks the correctness of the selected answer for the given chemistry question.
    It programmatically verifies the products of two Michael addition reactions based on IUPAC naming rules and reaction mechanisms.
    """
    
    # The multiple-choice options provided in the question
    options = {
        'A': {
            'A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'B': "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'B': {
            'A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'B': "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'C': {
            'A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'B': "3-methyl-4-nitrohexanenitrile"
        },
        'D': {
            'A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'B': "3-methyl-4-nitrohexanenitrile"
        }
    }

    # The final answer to be checked
    final_answer_choice = 'D'

    # --- Step 1: Determine the correct product for Reaction A ---
    # Reaction: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK)
    # Analysis: This is a Michael addition. The base (t-BuOK) deprotonates the ketone.
    # Constraint 1: The ketone has two alpha-carbons, C2 and C6. C2 is quaternary and has no alpha-protons.
    # C6 is tertiary and has one alpha-proton. Therefore, the enolate must form at C6.
    # Constraint 2: The C6 enolate attacks the beta-carbon of ethyl acrylate.
    # Constraint 3: Naming the product. The ester has higher priority than the ketone, so the parent chain is ethyl propanoate.
    # The cyclohexanone ring is a substituent at C3 of the propanoate.
    # Numbering the substituent ring starts at the point of attachment (original C6), which becomes C1'.
    # - C1' has a methyl group.
    # - C2' is the original carbonyl carbon (oxo).
    # - C3' is the original C2, with an ethyl and a methyl group.
    # This leads to the substituent name: (3-ethyl-1,3-dimethyl-2-oxocyclohexyl).
    correct_product_A = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"

    # --- Step 2: Determine the correct product for Reaction B ---
    # Reaction: 1-nitropropane + (KOH, (E)-but-2-enenitrile, H2O)
    # Analysis: This is a nitro-Michael addition.
    # Constraint 4: The base (KOH) deprotonates the carbon alpha to the nitro group in 1-nitropropane (a 3-carbon molecule).
    # Constraint 5: The resulting carbanion attacks the beta-carbon of (E)-but-2-enenitrile (a 4-carbon molecule).
    # Constraint 6: Naming the product. The principal functional group is the nitrile (-CN). The longest carbon chain including the nitrile is 6 carbons long (3 from nitropropane + 3 from butenenitrile, as one C is the nitrile C).
    # The parent name is hexanenitrile.
    # Numbering from the nitrile carbon as C1 gives: CN(1)-CH2(2)-CH(CH3)(3)-CH(NO2)(4)-CH2(5)-CH3(6).
    correct_product_B = "3-methyl-4-nitrohexanenitrile"

    # --- Step 3: Verify the selected answer ---
    if final_answer_choice not in options:
        return f"Error: The provided answer '{final_answer_choice}' is not a valid option (A, B, C, or D)."

    proposed_A = options[final_answer_choice]['A']
    proposed_B = options[final_answer_choice]['B']
    
    errors = []

    # Check if Product A from the selected option is correct
    if proposed_A != correct_product_A:
        error_msg = (
            f"Product A is incorrect. The proposed name is '{proposed_A}'.\n"
            f"Reason: The Michael addition on 2-ethyl-2,6-dimethylcyclohexan-1-one can only occur via the enolate at C6, as C2 has no alpha-protons. "
            f"The correct IUPAC name for the resulting product is '{correct_product_A}'."
        )
        errors.append(error_msg)

    # Check if Product B from the selected option is correct
    if proposed_B != correct_product_B:
        error_msg = (
            f"Product B is incorrect. The proposed name is '{proposed_B}'.\n"
            f"Reason: The Michael addition of a 3-carbon nucleophile (from 1-nitropropane) and a 4-carbon electrophile (but-2-enenitrile) results in a product with a 6-carbon parent chain. "
            f"The correct IUPAC name is '{correct_product_B}', a hexanenitrile derivative, not a butanenitrile."
        )
        errors.append(error_msg)

    if not errors:
        return "Correct"
    else:
        return "\n\n".join(errors)

# Execute the check and print the result
# Redirect stdout to capture the print output of the function
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

result = check_correctness()
print(result)

# Restore stdout
sys.stdout = old_stdout
output = captured_output.getvalue().strip()

# Final check to present the output clearly
if output == "Correct":
    print("Correct")
else:
    # This part is for displaying the reason for incorrectness
    print(output)
