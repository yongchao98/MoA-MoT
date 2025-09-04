def check_klinefelter_answer(selected_option: str):
    """
    Checks the correctness of the answer for the Klinefelter's vs. Down's syndrome question.

    The function verifies the answer against the core scientific principles:
    1. The mechanism must explain the consequence (phenotypic severity), not the cause (nondisjunction).
    2. The mechanism must be related to post-zygotic dosage compensation (X-inactivation).
    """

    # Define the properties of each option based on biological facts
    options_analysis = {
        'A': {
            'description': "chromatin methylation by histone methyltransferases in the post-zygote",
            'process_type': "Epigenetic Dosage Compensation (X-inactivation)",
            'timing': "Post-zygote",
            'relevance': "Explains consequence/severity",
            'is_correct': True
        },
        'B': {
            'description': "chiasmata resolution by separase in diakinesis",
            'process_type': "Meiosis",
            'timing': "Pre-zygote (Gametogenesis)",
            'relevance': "Explains cause of aneuploidy",
            'is_correct': False
        },
        'C': {
            'description': "progression of the polymerase alpha in the morula/blastocyst",
            'process_type': "DNA Replication",
            'timing': "Post-zygote",
            'relevance': "General cellular process, not specific to dosage compensation",
            'is_correct': False
        },
        'D': {
            'description': "attachment of spindle to kinetochores in the metaphase I",
            'process_type': "Meiosis",
            'timing': "Pre-zygote (Gametogenesis)",
            'relevance': "Explains cause of aneuploidy",
            'is_correct': False
        }
    }

    # Clean up the input answer
    answer_key = selected_option.strip().upper()

    if answer_key not in options_analysis:
        return f"Invalid option '{answer_key}'. Please provide one of {list(options_analysis.keys())}."

    # Check if the selected option meets the criteria for being correct
    if options_analysis[answer_key]['is_correct']:
        return "Correct"
    else:
        # Find the correct answer to provide a helpful error message
        correct_key = next(key for key, value in options_analysis.items() if value['is_correct'])
        
        # Formulate the reason for the error
        reason = (f"Incorrect. The selected answer '{answer_key}' is wrong.\n"
                  f"Reason: The question asks for a mechanism that explains the *less severe consequence* of Klinefelter's syndrome. "
                  f"Option {answer_key} describes '{options_analysis[answer_key]['process_type']}', which is related to the '{options_analysis[answer_key]['relevance']}'.\n"
                  f"The correct mechanism is X-chromosome inactivation, a form of dosage compensation that occurs post-zygotically. "
                  f"This corresponds to option '{correct_key}': '{options_analysis[correct_key]['description']}'.")
        return reason

# The final answer provided by the analysis is <<<A>>>.
final_answer_from_llm = "<<<A>>>"

# Extract the letter from the final answer format
answer_to_check = final_answer_from_llm.strip('<>').strip()

# Run the check
result = check_klinefelter_answer(answer_to_check)
print(result)