def check_chemistry_answer():
    """
    This function checks the correctness of the selected answer for the given
    organic chemistry question involving Michael additions.
    """

    # Define the correct products and reactant based on chemical principles.
    correct_A = "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate"
    correct_B = "3-(2-oxocyclohexyl)butanenitrile"
    correct_C = "cyclohexane-1,3-dione"

    # Define the options provided in the question.
    options = {
        "A": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "B": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        },
        "C": {
            "A": "trimethyl 2-(p-tolyl)propane-1,1,3-tricarboxylate",
            "B": "3-(2-oxocyclohexyl)butanenitrile",
            "C": "cyclohexane-1,3-dione"
        },
        "D": {
            "A": "trimethyl 3-(p-tolyl)propane-1,1,2-tricarboxylate",
            "B": "3-(2-hydroxycyclohex-1-en-1-yl)butanenitrile",
            "C": "2-hydroxycyclohexane-1,3-dione"
        }
    }

    # The answer provided by the LLM to be checked.
    llm_answer = "C"

    # Retrieve the components of the chosen answer.
    chosen_components = options.get(llm_answer)

    if not chosen_components:
        return f"Invalid answer key '{llm_answer}'. The key must be one of {list(options.keys())}."

    # Check each component against the correct answer.
    errors = []

    # Check Product A
    if chosen_components["A"] != correct_A:
        reason = (f"Product A is incorrect. The Michael addition of dimethyl malonate to methyl (E)-3-(p-tolyl)acrylate "
                  f"results in the nucleophile attacking the β-carbon, leading to '{correct_A}', not '{chosen_components['A']}'.")
        errors.append(reason)

    # Check Product B
    if chosen_components["B"] != correct_B:
        reason = (f"Product B is incorrect. The Stork enamine reaction followed by acidic workup yields the thermodynamically "
                  f"stable keto form, '{correct_B}', as the major product, not the enol tautomer '{chosen_components['B']}'.")
        errors.append(reason)

    # Check Reactant C
    if chosen_components["C"] != correct_C:
        reason = (f"Reactant C is incorrect. The Michael donor is '{correct_C}'. The name '{chosen_components['C']}' "
                  f"refers to the enol tautomer, not the standard name for the starting diketone reactant.")
        errors.append(reason)

    # Return the final verdict.
    if not errors:
        return "Correct"
    else:
        # Join all found errors into a single message.
        return "Incorrect. " + " ".join(errors)

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)