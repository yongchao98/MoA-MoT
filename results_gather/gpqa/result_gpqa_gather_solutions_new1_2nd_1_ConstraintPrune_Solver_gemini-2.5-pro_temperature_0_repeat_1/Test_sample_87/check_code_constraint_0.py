import re

def check_correctness():
    """
    This function checks the correctness of the provided answer to the immunology question.
    It models the key characteristics of the immunological processes listed in the options
    and compares them against the constraints described in the experimental scenario.
    """

    # The final answer provided by the LLM to be checked.
    llm_answer_text = "<<<D>>>"

    # Define the characteristics of each immunological process based on the options.
    # This represents the "knowledge base" for the checker.
    processes = {
        'A': {
            'name': 'class switching recombination',
            'location': 'secondary_lymphoid_organ',
            'timing': 'post_antigen_stimulation',
            'genetic_target': 'constant_region',
            'is_genetic_modification': True,
            'description': "Changes the antibody's function (isotype) by altering the constant region of the heavy chain."
        },
        'B': {
            'name': 'complement activation',
            'location': 'blood/tissue',
            'timing': 'post_antigen_stimulation',
            'genetic_target': 'none',
            'is_genetic_modification': False,
            'description': "A protein cascade in the blood that helps clear pathogens; does not involve genetic changes in B cells."
        },
        'C': {
            'name': 'VDJ recombination',
            'location': 'primary_lymphoid_organ',
            'timing': 'pre_antigen_encounter',
            'genetic_target': 'variable_region',
            'is_genetic_modification': True,
            'description': "Creates initial antibody diversity in developing B cells in the bone marrow, before antigen encounter."
        },
        'D': {
            'name': 'somatic hypermutation',
            'location': 'secondary_lymphoid_organ',
            'timing': 'post_antigen_stimulation',
            'genetic_target': 'variable_region',
            'is_genetic_modification': True,
            'description': "Introduces point mutations into the variable region of immunoglobulin genes in activated B cells to refine antibody affinity."
        }
    }

    # Define the constraints derived from the question's experimental scenario.
    question_constraints = {
        'location': {
            'value': 'secondary_lymphoid_organ',
            'reason': "The process occurs in Peyer's patches, which are secondary lymphoid organs."
        },
        'timing': {
            'value': 'post_antigen_stimulation',
            'reason': "It happens in a 'proliferating cell population' after antigen delivery, which is a post-stimulation event."
        },
        'genetic_target': {
            'value': 'variable_region',
            'reason': "The key finding is 'high variability' specifically in the 'variable heavy chain gene'."
        }
    }

    # Find the option that satisfies all constraints.
    correct_option_key = None
    for key, process_details in processes.items():
        is_match = True
        for constraint, details in question_constraints.items():
            if process_details.get(constraint) != details['value']:
                is_match = False
                break
        if is_match:
            correct_option_key = key
            break

    # Extract the letter from the LLM's answer string.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Could not parse the answer format: {llm_answer_text}"
    
    llm_answer_key = match.group(1)

    # Compare the LLM's answer with the derived correct answer.
    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness.
        llm_process = processes[llm_answer_key]
        correct_process = processes[correct_option_key]
        
        reason = (f"The provided answer '{llm_answer_key}' ({llm_process['name']}) is incorrect. "
                  f"The correct answer is '{correct_option_key}' ({correct_process['name']}).\n\n")
        
        reason += "Reasoning:\n"
        for constraint, details in question_constraints.items():
            correct_value = correct_process[constraint]
            llm_value = llm_process[constraint]
            if correct_value != llm_value:
                reason += f"- Constraint: {details['reason']}\n"
                reason += f"  - Correct process '{correct_process['name']}' satisfies this (value: '{correct_value}').\n"
                reason += f"  - Provided answer '{llm_process['name']}' fails this constraint (value: '{llm_value}').\n"
                reason += f"  - Specifically, {llm_process['name']} affects the '{llm_process['genetic_target']}' and/or occurs in a '{llm_process['location']}' at a time '{llm_process['timing']}'.\n"
        
        return reason

# Execute the check and print the result.
result = check_correctness()
print(result)