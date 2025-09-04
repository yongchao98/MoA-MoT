def check_immunology_question():
    """
    This function checks the correctness of the answer to the immunology question.
    It models the key features of the immunological processes and the conditions described in the question.
    """

    # Define the properties of each immunological process given in the options
    processes = {
        'A': {
            'name': 'Class Switching Recombination',
            'timing': 'During active immune response in mature B cells',
            'location': 'Secondary lymphoid organs (e.g., Peyer\'s patches)',
            'genetic_target': 'Constant region of the heavy chain',
            'effect': 'Changes antibody isotype (e.g., IgM to IgA), altering function but not antigen specificity.'
        },
        'B': {
            'name': 'Somatic Hypermutation (SHM)',
            'timing': 'During active immune response in mature B cells',
            'location': 'Secondary lymphoid organs (e.g., Peyer\'s patches)',
            'genetic_target': 'Variable region of heavy and light chains',
            'effect': 'Introduces point mutations, creating high variability to increase antigen affinity.'
        },
        'C': {
            'name': 'VDJ Recombination',
            'timing': 'During initial B cell development (before antigen encounter)',
            'location': 'Bone marrow',
            'genetic_target': 'Variable region of heavy and light chains',
            'effect': 'Creates the initial, diverse repertoire of B cell receptors.'
        },
        'D': {
            'name': 'Complement Activation',
            'timing': 'Innate immune response (can be activated by adaptive immunity)',
            'location': 'Blood and tissues',
            'genetic_target': 'None (it is a protein cascade, not a genetic modification process)',
            'effect': 'Opsonization, inflammation, and direct lysis of pathogens.'
        }
    }

    # Extract the key conditions from the question's scenario
    question_conditions = {
        'timing': 'During active immune response in mature B cells', # Implied by "proliferating cell population" after antigen challenge
        'location': 'Secondary lymphoid organs (e.g., Peyer\'s patches)', # Explicitly stated "Peyer patches"
        'genetic_target': 'Variable region of heavy and light chains', # Explicitly stated "variable heavy chain gene"
        'effect': 'Introduces point mutations, creating high variability to increase antigen affinity.' # Implied by "high variability" in a proliferating population
    }

    # The answer provided by the LLM
    llm_answer = 'B'

    # Check the provided answer against the conditions
    chosen_process = processes.get(llm_answer)

    if not chosen_process:
        return f"Invalid answer option '{llm_answer}'. Please choose from A, B, C, or D."

    # Perform the logical checks
    if chosen_process['timing'] != question_conditions['timing']:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because its timing is incorrect. "
                f"The process in the question occurs '{question_conditions['timing']}', but {chosen_process['name']} "
                f"occurs '{chosen_process['timing']}'.")

    if chosen_process['location'] != question_conditions['location']:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because its location is incorrect. "
                f"The process in the question occurs in '{question_conditions['location']}', but {chosen_process['name']} "
                f"occurs in the '{chosen_process['location']}'.")

    if chosen_process['genetic_target'] != question_conditions['genetic_target']:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because it targets the wrong gene region. "
                f"The question describes changes in the '{question_conditions['genetic_target']}', but {chosen_process['name']} "
                f"targets the '{chosen_process['genetic_target']}'.")
    
    # This is the most crucial check for this specific question
    if "high variability" not in chosen_process['effect'].lower():
         return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because its effect does not match the key observation. "
                f"The question observes 'high variability' in the variable gene, which is the hallmark of {processes['B']['name']}. "
                f"The effect of {chosen_process['name']} is '{chosen_process['effect']}'.")


    # If all conditions are met, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_immunology_question()
print(result)