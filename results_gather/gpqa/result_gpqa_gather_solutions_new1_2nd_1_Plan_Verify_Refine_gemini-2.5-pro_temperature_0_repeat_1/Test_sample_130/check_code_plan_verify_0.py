import re

def check_answer():
    """
    Checks the correctness of the LLM's answer by codifying the chemical reasoning.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_answer = "<<<D>>>"

    # --- Step 1: Define problem parameters and chemical knowledge ---

    # The options as provided in the question
    options = {
        'A': {'A 1H doublet at ~1.5 ppm', 'a 2H singlet at ~3.5 ppm'},
        'B': {'A 6H singlet at ~1 ppm', 'a 1H doublet at ~1.5 ppm'},
        'C': {'A 6H singlet at ~1 ppm', 'a 6H singlet at ~1.7 ppm'},
        'D': {'A 6H singlet at ~1.7 ppm', 'a 2H singlet at ~3.5 ppm'}
    }

    # Proton types and their expected NMR signal descriptions
    signal_assignments = {
        'H_anhydride': '2H singlet at ~3.5 ppm',
        'Me_vinyl': '6H singlet at ~1.7 ppm',
        'H_bridge': '1H doublet at ~1.5 ppm',
        'Me_bridgehead': '6H singlet at ~1.0 ppm'
    }

    # --- Step 2: Determine the major product based on stereoselectivity rules ---
    
    # The diene is 1,2,3,4-tetramethyl-1,3-cyclopentadiene, which is very bulky.
    # Rule: Extreme steric hindrance overrides the Alder-endo rule.
    # Therefore, the exo adduct is the major product.
    major_product_isomer = 'exo'
    
    reasoning_log = [
        "1. Stereoselectivity Rule: For the reaction between maleic anhydride and 1,2,3,4-tetramethyl-1,3-cyclopentadiene, the diene is exceptionally bulky.",
        "   - This steric hindrance overrides the standard Alder-endo rule.",
        "   - Conclusion: The major product is the 'exo' adduct."
    ]

    # --- Step 3: Determine the key NOESY interaction for the major product ---

    # In the 'exo' adduct, the anhydride protons (H_anhydride) are on the endo face,
    # spatially close to the vinylic methyl groups (Me_vinyl), which are also on the endo face.
    # In the 'endo' adduct, H_anhydride would be on the exo face, far from Me_vinyl.
    distinguishing_interaction_protons = ('H_anhydride', 'Me_vinyl')
    
    reasoning_log.append(
        "2. NOESY Analysis: A cross-peak present in the major product but absent in the minor indicates a unique spatial proximity."
    )
    reasoning_log.append(
        f"   - In the major ('{major_product_isomer}') product, the key spatial proximity is between {distinguishing_interaction_protons[0]} and {distinguishing_interaction_protons[1]}."
    )

    # --- Step 4: Map the interacting protons to their NMR signals ---
    
    signal1 = signal_assignments[distinguishing_interaction_protons[0]]
    signal2 = signal_assignments[distinguishing_interaction_protons[1]]
    expected_signals = {signal1, signal2}
    
    reasoning_log.append(
        f"3. Signal Assignment: The interacting protons correspond to the signals: '{signal1}' and '{signal2}'."
    )

    # --- Step 5: Find the matching option ---

    # Normalize strings for robust comparison (lowercase, remove articles)
    def normalize_signal(s):
        return re.sub(r'^(a|an)\s+', '', s.lower()).strip()

    normalized_expected_signals = {normalize_signal(s) for s in expected_signals}
    
    derived_correct_option = None
    for option_key, signal_set in options.items():
        normalized_option_signals = {normalize_signal(s) for s in signal_set}
        if normalized_option_signals == normalized_expected_signals:
            derived_correct_option = option_key
            break
            
    if not derived_correct_option:
        return "Error in checking logic: Could not match the derived interacting signals to any of the provided options."

    reasoning_log.append(
        f"4. Conclusion: The pair of signals matches option '{derived_correct_option}'."
    )

    # --- Step 6: Compare with the LLM's answer ---
    
    llm_option = llm_answer.strip('<>')
    
    if llm_option == derived_correct_option:
        return "Correct"
    else:
        error_message = (
            f"The provided answer '{llm_answer}' is incorrect.\n"
            f"The logically derived correct answer is '<<<{derived_correct_option}>>>'.\n\n"
            "--- REASONING ---\n"
        )
        error_message += "\n".join(reasoning_log)
        return error_message

# Execute the check
result = check_answer()
print(result)