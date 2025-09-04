import re

def check_chemistry_answer():
    """
    This function checks the correctness of the given answer by encoding the chemical properties
    of the compounds mentioned in the question.
    """

    # --- Chemical Facts Database ---
    # We encode the known properties of each compound based on chemical principles.
    compound_properties = {
        "benzoquinone": {
            "shows_tautomerism": False,
            "reason_tautomerism": "It lacks alpha-hydrogens, which are required for keto-enol tautomerism."
        },
        "cyclohexane-1,3,5-trione": {
            "shows_tautomerism": True,
            "reason_tautomerism": "It has alpha-hydrogens and can tautomerize into its stable aromatic enol form (phloroglucinol)."
        },
        "methyl 2-hydroxypropanoate": {
            "shows_optical_isomerism": True,
            "reason_optical": "Its second carbon is a chiral center, bonded to four different groups (H, OH, CH3, COOCH3)."
        },
        "dimethyl fumarate": {
            "shows_optical_isomerism": False,
            "reason_optical": "It is an achiral molecule with a plane of symmetry and no chiral centers."
        }
    }

    # --- Question Constraints ---
    # Part A asks for the compound that DOES NOT show tautomerism.
    # Part B asks for the compound that WILL show optical isomerism.
    
    # --- The LLM's Answer to be Checked ---
    # The provided answer is 'A', which corresponds to the text:
    # "A = benzoquinone, B = methyl 2-hydroxypropanoate"
    llm_answer_text = "A = benzoquinone, B = methyl 2-hydroxypropanoate"

    # --- Parsing the Answer ---
    try:
        # Use regex to safely extract compound names
        match = re.match(r"A = (.*), B = (.*)", llm_answer_text)
        if not match:
            return f"Invalid answer format. The answer '{llm_answer_text}' does not match the expected 'A = [compound], B = [compound]' format."
        
        compound_A_name = match.group(1).strip()
        compound_B_name = match.group(2).strip()

        # Verify that the parsed names are in our database
        if compound_A_name not in compound_properties or compound_B_name not in compound_properties:
             return f"Invalid compound name found in the answer: '{compound_A_name}' or '{compound_B_name}'."

    except Exception as e:
        return f"Failed to parse the answer string. Error: {e}"

    # --- Verification Logic ---
    # Check constraint for A: Must NOT show tautomerism.
    is_A_correct = compound_properties[compound_A_name]["shows_tautomerism"] is False
    if not is_A_correct:
        reason = compound_properties[compound_A_name]["reason_tautomerism"]
        return (f"Incorrect. The choice for A ('{compound_A_name}') is wrong. "
                f"The question asks for the compound that DOES NOT show tautomerism, but {compound_A_name} does. "
                f"Reason: {reason}")

    # Check constraint for B: Must show optical isomerism.
    is_B_correct = compound_properties[compound_B_name]["shows_optical_isomerism"] is True
    if not is_B_correct:
        reason = compound_properties[compound_B_name]["reason_optical"]
        return (f"Incorrect. The choice for B ('{compound_B_name}') is wrong. "
                f"The question asks for the compound that WILL show optical isomerism, but {compound_B_name} does not. "
                f"Reason: {reason}")

    # If both conditions are met
    return "Correct"

# Execute the check
result = check_chemistry_answer()
print(result)