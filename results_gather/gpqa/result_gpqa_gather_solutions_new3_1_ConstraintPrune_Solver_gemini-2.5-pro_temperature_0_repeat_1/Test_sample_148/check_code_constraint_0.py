import re

def check_answer_correctness(question_text: str, llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to a multiple-choice chemistry question.

    The function programmatically analyzes the constraints given in the question and
    evaluates whether the selected option satisfies all of them.

    Args:
        question_text: The string containing the question and options.
        llm_answer_text: The string containing the LLM's full response, including the final answer tag.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """
    # Extract the final answer choice (A, B, C, or D) from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: The answer format is incorrect. It should end with <<<A>>>, <<<B>>>, <<<C>>>, or <<<D>>>."
    
    llm_answer_choice = match.group(1)

    # --- Step 1: Define constraints based on the experimental data in the question ---
    # Constraint 1: From LC-MS showing two peaks with the same mass spectrum as the expected molecule.
    # This means the two species are isomers.
    is_isomeric = True
    
    # Constraint 2: From 1H NMR showing two distinct peaks for a single proton.
    # This means the species are distinguishable by standard NMR (i.e., not enantiomers).
    is_distinguishable_by_nmr = True
    
    # Constraint 3: From LC showing two clearly defined peaks.
    # This means the species are separable by standard chromatography (i.e., not enantiomers).
    is_separable_by_lc = True

    # --- Step 2: Define the properties of each possible answer option ---
    options_properties = {
        "A": {  # Diastereoisomers
            "name": "a mixture of diastereoisomers",
            "is_isomeric": True,
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True
        },
        "B": {  # Precursor contamination
            "name": "contamination with a precursor",
            "is_isomeric": False,  # A precursor has a different mass
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True
        },
        "C": {  # Enantiomers
            "name": "a mixture of enantiomers",
            "is_isomeric": True,
            "is_distinguishable_by_nmr": False,  # Indistinguishable in standard (achiral) NMR
            "is_separable_by_lc": False   # Not separable by standard (achiral) LC
        },
        "D": {  # 'Double coupling' product
            "name": "'double coupling'",
            "is_isomeric": False,  # A double-coupled product would have a higher mass
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True
        }
    }

    # --- Step 3: Determine the correct answer by checking which option satisfies all constraints ---
    correct_choice = None
    for choice, properties in options_properties.items():
        if (properties["is_isomeric"] == is_isomeric and
            properties["is_distinguishable_by_nmr"] == is_distinguishable_by_nmr and
            properties["is_separable_by_lc"] == is_separable_by_lc):
            correct_choice = choice
            break

    # --- Step 4: Compare the LLM's answer with the determined correct answer ---
    if llm_answer_choice == correct_choice:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrectness
        llm_props = options_properties[llm_answer_choice]
        reason = f"The answer '{llm_answer_choice}' ({llm_props['name']}) is incorrect.\n"
        
        # Check which specific constraints are violated
        if llm_props["is_isomeric"] != is_isomeric:
            reason += "- Constraint violated: The species must be isomers (same mass). The MS data shows both peaks have the same mass as the expected product, but your answer implies a different mass.\n"
        if llm_props["is_distinguishable_by_nmr"] != is_distinguishable_by_nmr:
            reason += "- Constraint violated: The species must be distinguishable by NMR. The NMR data shows two distinct peaks, but your answer implies they would give a single set of peaks in a standard experiment.\n"
        if llm_props["is_separable_by_lc"] != is_separable_by_lc:
            reason += "- Constraint violated: The species must be separable by LC. The LC data shows two distinct peaks, but your answer implies they would not separate on a standard column.\n"
            
        reason += f"The correct answer is '{correct_choice}' ({options_properties[correct_choice]['name']}), as it is the only option consistent with all experimental data (isomers that are separable by LC and distinguishable by NMR)."
        return reason

# You can run the check by passing the question and the LLM's answer to the function.
# For example:
# question = "..." (The full question text)
# llm_answer = "..." (The full answer text provided by the LLM)
# result = check_answer_correctness(question, llm_answer)
# print(result)