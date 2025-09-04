def check_answer_correctness():
    """
    This function checks the correctness of the answer to the immunology question.
    It models the key facts from the question and the properties of the biological processes.
    """

    # Define the key characteristics of each process listed in the options.
    # This represents the ground truth from immunology textbooks.
    processes = {
        "A": {
            "name": "VDJ recombination",
            "timing": "pre-antigen",  # Occurs during B cell development
            "location": "primary_lymphoid_organ",  # e.g., Bone Marrow
            "gene_region_affected": "variable",
            "purpose": "Creates initial receptor diversity"
        },
        "B": {
            "name": "complement activation",
            "timing": "post-antigen",
            "location": "blood/tissues",
            "gene_region_affected": "none",  # It's a protein cascade, not a genetic modification of lymphocytes
            "purpose": "Helps clear pathogens"
        },
        "C": {
            "name": "class switching recombination",
            "timing": "post-antigen",
            "location": "secondary_lymphoid_organ",  # e.g., Germinal centers in Peyer's patches
            "gene_region_affected": "constant",  # This is a key differentiator
            "purpose": "Changes antibody effector function (isotype)"
        },
        "D": {
            "name": "somatic hypermutation",
            "timing": "post-antigen",
            "location": "secondary_lymphoid_organ",  # e.g., Germinal centers in Peyer's patches
            "gene_region_affected": "variable",  # This is the key observation
            "purpose": "Introduces point mutations to increase antibody affinity"
        }
    }

    # Define the facts and constraints derived from the question's text.
    question_constraints = {
        "timing": "post-antigen",  # The response happens after antigen delivery.
        "location": "secondary_lymphoid_organ",  # The response is in Peyer's patches.
        "gene_region_affected": "variable",  # The key finding is "high variability in the variable heavy chain gene".
    }

    # The final answer provided by the LLM being checked.
    llm_answer = "D"

    # Determine the correct answer by finding which process matches all constraints.
    correct_option = None
    for option, properties in processes.items():
        if all(properties[key] == question_constraints[key] for key in question_constraints):
            correct_option = option
            break
    
    # If the LLM's answer matches the derived correct answer, it is correct.
    if llm_answer == correct_option:
        return "Correct"
    else:
        # If not, generate a reason for why the LLM's answer is wrong.
        wrong_answer_properties = processes.get(llm_answer)
        if not wrong_answer_properties:
            return f"Invalid answer: The option '{llm_answer}' was not one of the choices."

        reasons = []
        if wrong_answer_properties["timing"] != question_constraints["timing"]:
            reasons.append(f"the process occurs {wrong_answer_properties['timing'].replace('_', ' ')}, but the experiment is {question_constraints['timing'].replace('_', ' ')}")
        
        if wrong_answer_properties["location"] != question_constraints["location"]:
            reasons.append(f"the process occurs in a {wrong_answer_properties['location'].replace('_', ' ')}, but the experiment takes place in a {question_constraints['location'].replace('_', ' ')} (Peyer's patch)")

        if wrong_answer_properties["gene_region_affected"] != question_constraints["gene_region_affected"]:
            reasons.append(f"the process affects the '{wrong_answer_properties['gene_region_affected']}' gene region, but the question specifies variability in the '{question_constraints['gene_region_affected']}' region")

        return f"Incorrect. The answer '{llm_answer}' ({wrong_answer_properties['name']}) is wrong because " + ", and ".join(reasons) + "."

# Execute the check and print the result.
result = check_answer_correctness()
print(result)