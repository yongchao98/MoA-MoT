import re

def check_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the Diels-Alder NOESY problem.

    The logic is as follows:
    1.  **Identify Reactants and Reaction:** The reaction is a Diels-Alder cycloaddition between maleic anhydride and 1,2,3,4-tetramethyl-1,3-cyclopentadiene.
    2.  **Identify Major Product:** According to the Alder endo rule, the major product formed under kinetic control is the 'endo' adduct.
    3.  **Identify Key Spatial Proximity for NOESY:** The question asks for a cross-peak present in the major product but absent in the minor. This requires finding a pair of protons that are close in the 'endo' adduct but far apart in the 'exo' adduct.
    4.  **Analyze Stereochemistry:** In the 'endo' adduct, the anhydride ring is tucked under the C=C double bond. This brings the anhydride protons very close to the vinylic methyl groups. This proximity is absent in the 'exo' adduct, where the anhydride ring points away from the double bond.
    5.  **Assign NMR Signals:** Based on the options and standard chemical shifts:
        - Anhydride protons: "a 2H singlet at ~3.5 ppm"
        - Vinylic methyl protons: "A 6H singlet at ~1.7 ppm"
        - Bridgehead methyl protons: "A 6H singlet at ~1 ppm"
        - Bridge proton: "a 1H doublet at ~1.5 ppm"
    6.  **Conclude the Correct Interaction:** The distinguishing cross-peak must be between the anhydride protons and the vinylic methyl protons.
    7.  **Match to Options:** This interaction corresponds to the pair of signals: "a 2H singlet at ~3.5 ppm" and "A 6H singlet at ~1.7 ppm". This is option C.
    8.  **Compare:** The derived correct option is compared with the LLM's answer.
    """

    # Step 1: Define proton signals based on the options provided in the question.
    # We use labels for clarity and map them to the exact text in the options.
    proton_signals = {
        "anhydride_H": "a 2H singlet at ~3.5 ppm",
        "vinylic_Me": "A 6H singlet at ~1.7 ppm",
        "bridgehead_Me": "A 6H singlet at ~1 ppm",
        "bridge_H": "a 1H doublet at ~1.5 ppm"
    }

    # Step 2: Determine the correct interaction based on chemical principles.
    # The major product is the 'endo' adduct. The unique NOESY cross-peak for the
    # endo adduct is between the anhydride protons and the vinylic methyl groups.
    correct_interaction_labels = {"anhydride_H", "vinylic_Me"}
    correct_interaction_descriptions = {proton_signals[label] for label in correct_interaction_labels}

    # Step 3: Define the options from the question as sets of descriptions for easy comparison.
    options = {
        "A": {"A 6H singlet at ~1 ppm", "a 1H doublet at ~1.5 ppm"},
        "B": {"A 6H singlet at ~1 ppm", "a 6H singlet at ~1.7 ppm"},
        "C": {"A 6H singlet at ~1.7 ppm", "a 2H singlet at ~3.5 ppm"},
        "D": {"a 1H doublet at ~1.5 ppm", "a 2H singlet at ~3.5 ppm"}
    }

    # Step 4: Find the letter of the correct option by matching the set of descriptions.
    correct_option_letter = None
    for letter, descriptions in options.items():
        if descriptions == correct_interaction_descriptions:
            correct_option_letter = letter
            break

    if correct_option_letter is None:
        # This is a sanity check for the checker's logic.
        return "Internal checker error: Could not map the correct chemical interaction to any of the options."

    # Step 5: Extract the LLM's chosen answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The provided answer text does not contain a final answer in the format <<<X>>>."

    llm_choice = match.group(1)

    # Step 6: Compare the LLM's choice with the correct answer and return the result.
    if llm_choice == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{llm_choice}' is incorrect. The correct answer is '{correct_option_letter}'.\n"
            "Reasoning:\n"
            "1. The major product of the Diels-Alder reaction is the 'endo' adduct.\n"
            "2. In the 'endo' adduct, the anhydride protons are spatially close to the vinylic methyl groups, while they are far apart in the minor 'exo' adduct.\n"
            "3. This unique proximity in the major product gives rise to the distinguishing NOESY cross-peak.\n"
            f"4. The signals for these protons are '{proton_signals['anhydride_H']}' and '{proton_signals['vinylic_Me']}'.\n"
            f"5. This pair of signals corresponds to option {correct_option_letter}."
        )
        return reason

# Example of how to use the function with the provided answer text:
# llm_answer = """
# ... (paste the full text of the answer to be checked here) ...
# <<<C>>>
# """
# result = check_correctness(llm_answer)
# print(result)