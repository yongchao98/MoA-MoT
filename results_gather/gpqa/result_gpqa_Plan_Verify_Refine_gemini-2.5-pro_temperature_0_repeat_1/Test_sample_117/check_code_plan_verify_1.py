def check_chemistry_answer():
    """
    This function verifies the correctness of an LLM's answer to a specific
    organic chemistry question by encoding the reaction rules.

    Question: What is the major outcome of the reaction between 4,4-dimethylcyclopent-1-enol and bromine?
    LLM's Answer Choice: B
    LLM's Reasoning: The reaction of an enol with bromine results in an α-bromo ketone...
                     The bromine will add to the α-carbon (C2), and the enol will
                     tautomerize to a ketone, placing a carbonyl group at C1.
                     The resulting product is 2-bromo-4,4-dimethylcyclopentanone.
    """

    # --- Define the problem parameters and the LLM's answer ---
    reactant = "4,4-dimethylcyclopent-1-enol"
    reagent = "bromine"
    llm_answer_choice = "B"
    llm_reasoning = """
    The reaction of an enol with bromine results in an α-bromo ketone. For 4,4-dimethylcyclopent-1-enol,
    the double bond is between C1 (with the -OH) and C2. The bromine will add to the α-carbon (C2),
    and the enol will tautomerize to a ketone, placing a carbonyl group at C1. The resulting product is
    2-bromo-4,4-dimethylcyclopentanone.
    """
    options = {
        "A": "4-bromo-4,4-dimethylcyclopentanone",
        "B": "2-bromo-4,4-dimethylcyclopentanone",
        "C": "(1R,2R)-1,2-dibromo-4,4-dimethylcyclopentanol",
        "D": "(1R,2S)-1,2-dibromo-4,4-dimethylcyclopentanol"
    }
    llm_product_name = options.get(llm_answer_choice)

    # --- Apply Chemical Rules to Predict the Correct Product ---

    # Rule 1: The reaction is the bromination of an enol.
    # Rule 2: The product is an alpha-bromo ketone.
    # Analysis of reactant '4,4-dimethylcyclopent-1-enol':
    # - The enol group involves C1 (with -OH) and C2.
    # - The alpha-carbon relative to the hydroxyl is C2.
    # Prediction:
    # - The enol's -OH at C1 becomes a ketone (C=O). The parent name becomes 'cyclopentanone'.
    # - Bromine adds to the alpha-carbon, C2. This gives a '2-bromo' substituent.
    # - The '4,4-dimethyl' group is unaffected.
    # - Assembling the name alphabetically: 2-bromo-4,4-dimethylcyclopentanone.
    expected_product_name = "2-bromo-4,4-dimethylcyclopentanone"

    # --- Verification Step ---

    # 1. Check if the chosen product name is correct.
    if llm_product_name != expected_product_name:
        return (f"Incorrect. The chosen product '{llm_product_name}' is wrong. "
                f"The correct product from the alpha-halogenation of the enol is '{expected_product_name}'.")

    # 2. Check if the reasoning is sound and covers the key steps.
    if "α-bromo ketone" not in llm_reasoning and "alpha-bromo ketone" not in llm_reasoning:
        return "Incorrect. The reasoning fails to identify the correct product class, which is an α-bromo ketone."
    
    if "add to the α-carbon (C2)" not in llm_reasoning:
        return "Incorrect. The reasoning fails to specify that bromine adds to the alpha-carbon (C2)."

    if "tautomerize to a ketone" not in llm_reasoning or "carbonyl group at C1" not in llm_reasoning:
        return "Incorrect. The reasoning is incomplete as it does not fully explain the tautomerization to a ketone at C1."

    # If all checks pass, the answer is correct.
    return "Correct"

# To check the answer, you would run this function.
# The function will return "Correct" if the LLM's answer and reasoning are valid.
# print(check_chemistry_answer())