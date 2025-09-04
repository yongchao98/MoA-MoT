def check_answer():
    """
    This function checks the correctness of the provided LLM's answer
    regarding the optical activity of a list of chemical compounds.
    """
    # Ground truth analysis based on chemical principles.
    correct_analysis = {
        "(Z)-1-chloro-2-methylbut-1-ene": {
            "is_active": False,
            "reason": "Achiral due to a plane of symmetry containing the C=C double bond."
        },
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": {
            "is_active": True,
            "reason": "Chiral. The name specifies a single enantiomer."
        },
        "(2R,3S)-2,3-dimethylsuccinic acid": {
            "is_active": False,
            "reason": "Achiral. This is a meso compound with an internal plane of symmetry."
        },
        "(2R,3R)-2,3-dimethylsuccinic acid": {
            "is_active": True,
            "reason": "Chiral. This is a specific enantiomer of a chiral molecule."
        },
        "(R)-cyclohex-3-en-1-ol": {
            "is_active": True,
            "reason": "Chiral. C1 is a stereocenter and a single (R) enantiomer is specified."
        },
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": {
            "is_active": False,
            "reason": "Achiral due to high symmetry (C3 axis and planes of symmetry)."
        },
        "1-cyclopentyl-3-methylbutan-1-one": {
            "is_active": False,
            "reason": "Achiral. Lacks a stereocenter and has a plane of symmetry."
        }
    }

    # Extract the analysis from the LLM's provided answer.
    llm_analysis = {
        "(Z)-1-chloro-2-methylbut-1-ene": False,
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": True,
        "(2R,3S)-2,3-dimethylsuccinic acid": False,
        "(2R,3R)-2,3-dimethylsuccinic acid": True,
        "(R)-cyclohex-3-en-1-ol": True,
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": False,
        "1-cyclopentyl-3-methylbutan-1-one": False
    }
    llm_final_answer_option = "D"

    # --- Verification Step ---
    errors = []
    
    # 1. Verify individual compound assessments
    for name, data in correct_analysis.items():
        if name not in llm_analysis:
            errors.append(f"Error: The LLM's answer did not include an analysis for '{name}'.")
            continue
        
        if llm_analysis[name] != data["is_active"]:
            expected = "active" if data["is_active"] else "inactive"
            got = "active" if llm_analysis[name] else "inactive"
            errors.append(f"Error on '{name}': Expected optical activity to be '{expected}', but the LLM stated it is '{got}'. Reason: {data['reason']}")

    # 2. Verify the final count
    correct_count = sum(1 for data in correct_analysis.values() if data["is_active"])
    llm_count = sum(1 for active in llm_analysis.values() if active)

    if correct_count != llm_count:
        errors.append(f"Error in final count: The correct number of optically active compounds is {correct_count}, but the LLM's answer implies a count of {llm_count}.")

    # 3. Verify the multiple-choice option
    answer_map = {'A': 4, 'B': 5, 'C': 2, 'D': 3}
    if llm_final_answer_option not in answer_map:
        errors.append(f"Error: The final answer option '{llm_final_answer_option}' is not a valid choice.")
    elif answer_map[llm_final_answer_option] != correct_count:
        errors.append(f"Error in final answer choice: The option '{llm_final_answer_option}' corresponds to a count of {answer_map[llm_final_answer_option]}, but the correct count is {correct_count}.")

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check
result = check_answer()
print(result)