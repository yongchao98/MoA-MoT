import re

def check_chemistry_answer(llm_answer):
    """
    Checks the correctness of the LLM's answer to the organic chemistry question.

    The function breaks down the problem into key logical steps:
    1.  Reactant Identification: Confirms the reactants are maleic anhydride and a substituted cyclopentadiene.
    2.  Stereochemical Prediction: Evaluates the most critical point - whether the major product is the 'endo' or 'exo' isomer. It correctly identifies that for this highly substituted diene, steric hindrance overrides the standard 'endo rule', making the 'exo' isomer the major product.
    3.  NOESY Analysis: Determines the unique spatial proximities for both isomers.
        - In the 'exo' isomer, the anhydride protons are close to the vinylic methyl protons.
        - In the 'endo' isomer, the anhydride protons are close to the C7 bridge protons.
    4.  Conclusion: Combines the prediction (exo is major) with the NOESY analysis to determine the expected cross-peak in the major product's spectrum.
    5.  Verification: Compares this expected outcome with the provided answer option.
    """
    
    # Step 1: Define the problem's components and the final answer to check
    try:
        final_answer = llm_answer.split('<<<')[-1].split('>>>')[0].strip()
    except IndexError:
        return "Invalid answer format. The answer should be enclosed in '<<< >>>'."

    # Define the proton signals and their assignments
    signals = {
        "anhydride_H": "2H singlet at ~3.5 ppm",
        "vinylic_Me": "6H singlet at ~1.7 ppm",
        "bridgehead_Me": "6H singlet at ~1.0 ppm",
        "bridge_H": "1H doublet at ~1.5 ppm"
    }

    # Define the options and the interactions they represent
    options = {
        "A": (signals["bridge_H"], signals["anhydride_H"]),
        "B": (signals["bridgehead_Me"], signals["vinylic_Me"]),
        "C": (signals["bridgehead_Me"], signals["bridge_H"]),
        "D": (signals["vinylic_Me"], signals["anhydride_H"])
    }

    # Step 2: Determine the major product (Exo vs. Endo)
    # The diene is 1,2,3,4-tetramethyl-1,3-cyclopentadiene. This is extremely bulky.
    # The standard "endo rule" is based on electronic effects (secondary orbital overlap).
    # However, severe steric hindrance can override this rule.
    # The endo approach would cause a major steric clash between the anhydride and the four methyl groups.
    # The exo approach is sterically much more favorable.
    # Conclusion: The major product is the EXO adduct.
    major_product_isomer = "exo"
    
    # Step 3: Determine the key NOESY interaction for each isomer
    # In the EXO isomer, the anhydride ring is on the same face as the vinylic methyl groups.
    exo_interaction = (signals["vinylic_Me"], signals["anhydride_H"])
    
    # In the ENDO isomer, the anhydride ring is on the same face as the C7 bridge.
    endo_interaction = (signals["bridge_H"], signals["anhydride_H"])

    # Step 4: Identify the interaction expected for the major product
    if major_product_isomer == "exo":
        expected_interaction = exo_interaction
        reasoning = "The major product is the 'exo' isomer due to overwhelming steric hindrance from the four methyl groups on the diene. In the 'exo' isomer, the anhydride protons (~3.5 ppm) are spatially close to the vinylic methyl protons (~1.7 ppm)."
    else: # major_product_isomer == "endo"
        expected_interaction = endo_interaction
        reasoning = "The major product is the 'endo' isomer based on the standard endo rule. In the 'endo' isomer, the anhydride protons (~3.5 ppm) are spatially close to the C7 bridge protons (~1.5 ppm)."

    # Step 5: Check if the provided answer matches the expected interaction
    if final_answer not in options:
        return f"Invalid option '{final_answer}'. The answer must be one of {list(options.keys())}."

    # The order of signals in the tuple doesn't matter for the interaction
    provided_interaction = options[final_answer]
    is_correct = (set(provided_interaction) == set(expected_interaction))

    if is_correct:
        return "Correct"
    else:
        # Find which option *should* have been correct
        correct_option = None
        for opt, interaction in options.items():
            if set(interaction) == set(expected_interaction):
                correct_option = opt
                break
        
        error_message = f"Incorrect. The provided answer is {final_answer}, but the correct answer is {correct_option}.\n"
        error_message += f"Reasoning: {reasoning}\n"
        error_message += f"The question asks for the NOESY cross-peak in the major product. {reasoning} This corresponds to the interaction between the '{expected_interaction[0]}' and the '{expected_interaction[1]}', which is option {correct_option}."
        
        return error_message

# The final answer from the LLM response
llm_final_answer = "<<<D>>>"

# Run the check
result = check_chemistry_answer(llm_final_answer)
print(result)