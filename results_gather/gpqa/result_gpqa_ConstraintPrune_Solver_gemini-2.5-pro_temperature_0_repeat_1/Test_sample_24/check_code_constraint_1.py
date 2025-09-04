def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for a chemistry question
    by verifying the functional groups of reactants for two name reactions.
    """
    
    # The answer provided by the LLM to be checked
    llm_answer = "C"

    # The options given in the question
    options = {
        'A': {
            'A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'B': '4-methyl-1-phenylpent-3-en-1-one'
        },
        'B': {
            'A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        },
        'C': {
            'A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        },
        'D': {
            'A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'B': '4-methyl-1-phenylpent-3-en-1-one'
        }
    }

    def get_functional_group(name):
        """Identifies the primary functional group from the IUPAC name based on suffixes/keywords."""
        name_lower = name.lower()
        # Check for 'diol' first as it's more specific than 'ol'
        if 'diol' in name_lower:
            return 'diol'
        if 'ol' in name_lower:
            return 'alcohol'
        if 'one' in name_lower:
            return 'ketone'
        # 'oxy' indicates an ether linkage (R-O-R')
        if 'oxy' in name_lower:
            return 'ether'
        return 'unknown'

    # Define the constraints for each reaction
    # Reaction A (Pinacol rearrangement) requires a diol.
    # Reaction B ([1,2]-Wittig rearrangement) requires an ether.
    
    valid_options = []
    analysis_log = {}

    for option_key, reactants in options.items():
        reactant_A_name = reactants['A']
        reactant_B_name = reactants['B']

        fg_A = get_functional_group(reactant_A_name)
        is_A_valid = (fg_A == 'diol')
        
        fg_B = get_functional_group(reactant_B_name)
        is_B_valid = (fg_B == 'ether')

        analysis_log[option_key] = {
            'A_valid': is_A_valid,
            'A_fg': fg_A,
            'B_valid': is_B_valid,
            'B_fg': fg_B
        }

        if is_A_valid and is_B_valid:
            valid_options.append(option_key)

    # Evaluate the LLM's answer
    if llm_answer not in options:
        return f"Incorrect. The provided answer '{llm_answer}' is not one of the possible options: {list(options.keys())}."

    is_llm_answer_correct = (len(valid_options) == 1 and valid_options[0] == llm_answer)

    if is_llm_answer_correct:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness
        llm_option_analysis = analysis_log[llm_answer]
        reasons = []
        if not llm_option_analysis['A_valid']:
            reasons.append(f"Reaction A requires a 'diol' for the Pinacol rearrangement, but in option {llm_answer}, reactant A is a(n) '{llm_option_analysis['A_fg']}'.")
        if not llm_option_analysis['B_valid']:
            reasons.append(f"Reaction B requires an 'ether' for the [1,2]-Wittig rearrangement, but in option {llm_answer}, reactant B is a(n) '{llm_option_analysis['B_fg']}'.")
        
        if len(valid_options) == 0:
            return f"Incorrect. The provided answer '{llm_answer}' is wrong, and in fact, no option satisfies both conditions. For option {llm_answer}: " + " ".join(reasons)
        else:
            correct_answer = valid_options[0]
            return f"Incorrect. The correct answer is {correct_answer}, not {llm_answer}. The provided answer fails because: " + " ".join(reasons)

# Execute the check and print the result
result = check_chemistry_answer()
print(result)