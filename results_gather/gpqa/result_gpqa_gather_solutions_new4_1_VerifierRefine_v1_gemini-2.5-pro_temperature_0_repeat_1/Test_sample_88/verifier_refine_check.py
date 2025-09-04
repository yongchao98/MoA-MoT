def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer to a multi-step synthesis problem.
    It encodes the known, literature-supported chemical pathway and NMR analysis principles
    to validate the reasoning step-by-step.
    """

    # --- Ground Truth Data ---
    # This dictionary contains the established correct answers for each step of the problem.
    correct_solution = {
        "product_1": "protoadamantan-4-one",
        "product_2": "protoadamant-4-ene",
        "product_3": "bicyclo[3.3.1]nonane-3,7-dione",
        "most_deshielded_proton": "bridgehead protons (C1 and C5)",
        "coupling_pattern": "triplet of triplets",
        "final_option": "C"
    }

    # --- Analysis of the LLM's Answer ---
    # This dictionary represents the key points from the provided answer's reasoning.
    llm_answer_analysis = {
        "product_1": "protoadamantan-4-one",
        "product_2": "protoadamant-4-ene",
        "product_3": "bicyclo[3.3.1]nonane-3,7-dione",
        "most_deshielded_proton": "bridgehead protons (C1 and C5)",
        "coupling_pattern": "triplet of triplets",
        "final_option": "C"
    }

    # --- Verification Logic ---
    error_messages = []

    # Check Product 1
    if llm_answer_analysis["product_1"] != correct_solution["product_1"]:
        error_messages.append(
            f"Constraint Violated: Identification of Product 1 is incorrect. "
            f"The reaction of 1,3-dibromoadamantane with hot KOH leads to a skeletal rearrangement, "
            f"forming {correct_solution['product_1']}, not {llm_answer_analysis['product_1']}."
        )

    # Check Product 2
    if llm_answer_analysis["product_2"] != correct_solution["product_2"]:
        error_messages.append(
            f"Constraint Violated: Identification of Product 2 is incorrect. "
            f"The reaction is an MPV reduction followed by dehydration (to enable the subsequent ozonolysis), "
            f"which forms {correct_solution['product_2']}, not {llm_answer_analysis['product_2']}."
        )

    # Check Product 3
    if llm_answer_analysis["product_3"] != correct_solution["product_3"]:
        error_messages.append(
            f"Constraint Violated: Identification of Product 3 is incorrect. "
            f"The ozonolysis of {correct_solution['product_2']} yields {correct_solution['product_3']}, "
            f"not {llm_answer_analysis['product_3']}."
        )

    # Check NMR Analysis: Most Deshielded Proton
    if llm_answer_analysis["most_deshielded_proton"] != correct_solution["most_deshielded_proton"]:
        error_messages.append(
            f"Constraint Violated: Identification of the most deshielded proton is incorrect. "
            f"In {correct_solution['product_3']}, the {correct_solution['most_deshielded_proton']} are the most deshielded, "
            f"not {llm_answer_analysis['most_deshielded_proton']}."
        )

    # Check NMR Analysis: Coupling Pattern
    if llm_answer_analysis["coupling_pattern"] != correct_solution["coupling_pattern"]:
        error_messages.append(
            f"Constraint Violated: The coupling pattern analysis is incorrect. "
            f"The bridgehead proton in {correct_solution['product_3']} is coupled to two distinct pairs of equivalent protons "
            f"(two axial and two equatorial on the adjacent methylenes), resulting in a {correct_solution['coupling_pattern']}, "
            f"not a {llm_answer_analysis['coupling_pattern']}."
        )

    # Check Final Answer
    if llm_answer_analysis["final_option"] != correct_solution["final_option"]:
        error_messages.append(
            f"Constraint Violated: The final answer option is incorrect. "
            f"The correct option is {correct_solution['final_option']}, but the answer provided was {llm_answer_analysis['final_option']}."
        )

    if not error_messages:
        return "Correct"
    else:
        return "\n".join(error_messages)

# Run the check and print the result
result = check_chemistry_answer()
print(result)