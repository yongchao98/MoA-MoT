import sys
from io import StringIO

def check_chemistry_answer():
    """
    This function checks the correctness of the selected answer for a chemistry question involving two Michael addition reactions.
    It determines the correct IUPAC names for both products based on chemical principles and compares them to the names in the chosen option.
    """
    
    # Define the options provided in the question
    options = {
        'A': {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate", 
            "B": "3-methyl-4-nitrohexanenitrile"
        },
        'B': {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate", 
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'C': {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate", 
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'D': {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate", 
            "B": "3-methyl-4-nitrohexanenitrile"
        }
    }

    # The answer to be checked is 'A'
    selected_option_letter = 'A'

    # --- Step 1: Determine the correct product name for Reaction A ---
    # Reaction: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK)
    # 1. Identify alpha-protons: The ketone has alpha-carbons at C2 and C6.
    #    - C2 is quaternary and has no alpha-protons.
    #    - C6 is tertiary and has one alpha-proton.
    # 2. Enolate formation: The base (t-BuOK) can only deprotonate C6.
    # 3. Michael addition: The C6 enolate attacks ethyl acrylate.
    # 4. IUPAC Naming: The ester has higher priority than the ketone.
    #    - Parent chain: ethyl propanoate.
    #    - The cyclohexanone ring is a substituent at C3 of the propanoate.
    #    - Numbering the substituent ring: The attachment point (original C6) is C1'. The carbonyl (original C1) is C2'.
    #    - Substituents on the ring: methyl at C1', oxo at C2', ethyl and methyl at C3'.
    #    - Full name: (3-ethyl-1,3-dimethyl-2-oxocyclohexyl)
    correct_name_A = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"

    # --- Step 2: Determine the correct product name for Reaction B ---
    # Reaction: 1-nitropropane + (E)-but-2-enenitrile (KOH)
    # 1. Nucleophile: The alpha-carbon of 1-nitropropane is deprotonated by KOH.
    # 2. Michael addition: The resulting carbanion attacks (E)-but-2-enenitrile.
    # 3. Product structure: CH3-CH2-CH(NO2)-CH(CH3)-CH2-CN
    # 4. IUPAC Naming: The nitrile is the principal functional group.
    #    - Parent chain: The longest chain including the nitrile carbon has 6 carbons -> hexanenitrile.
    #    - Numbering from CN as C1: CN(1)-CH2(2)-CH(CH3)(3)-CH(NO2)(4)-CH2(5)-CH3(6).
    #    - Substituents: methyl at C3, nitro at C4.
    correct_name_B = "3-methyl-4-nitrohexanenitrile"

    # --- Step 3: Validate the selected answer ---
    errors = []
    
    # Check if the selected option exists
    if selected_option_letter not in options:
        return f"Invalid option '{selected_option_letter}' provided. Valid options are {list(options.keys())}."

    proposed_answer = options[selected_option_letter]
    proposed_name_A = proposed_answer.get("A")
    proposed_name_B = proposed_answer.get("B")

    # Check Product A's name
    if proposed_name_A != correct_name_A:
        reason = (f"Product A is incorrect. The correct name is '{correct_name_A}'. "
                  f"The proposed name '{proposed_name_A}' is wrong. The enolate must form at C6, "
                  "and IUPAC naming of the resulting substituent places the oxo group at position 2, not 4 or 5.")
        errors.append(reason)

    # Check Product B's name
    if proposed_name_B != correct_name_B:
        reason = (f"Product B is incorrect. The correct name is '{correct_name_B}'. "
                  f"The proposed name '{proposed_name_B}' is wrong. The Michael addition of a 3-carbon unit (from 1-nitropropane) "
                  "to a 4-carbon unit (but-2-enenitrile) results in a product with a 6-carbon parent chain (hexanenitrile), "
                  "not a 4-carbon chain (butanenitrile).")
        errors.append(reason)

    # --- Step 4: Return the result ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. " + " ".join(errors)

# Execute the check and print the result
result = check_chemistry_answer()
print(result)