import re

def check_answer():
    """
    Checks the correctness of the LLM's answer for the 1H NMR problem.
    """
    # --- Data from the question ---
    nmr_data = {
        'signal_1': {'ppm': 7.0, 'integration': 1, 'multiplicity': 'd', 'J_Hz': 16.0},
        'signal_2': {'ppm': 5.5, 'integration': 1, 'multiplicity': 'dq'},
        'signal_3': {'ppm': 2.1, 'integration': 3, 'multiplicity': 's'},
        'signal_4': {'ppm': 1.6, 'integration': 3, 'multiplicity': 'd'},
    }

    options = {
        'A': 'Trans-butenyl acetate',
        'B': 'Trans-propenyl acetate',
        'C': 'Cis-propenyl acetate',
        'D': 'Cis-butenyl acetate'
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_text = "<<<B>>>"
    
    # --- Analysis Logic ---
    
    # Step 1: Check the total proton count to determine the carbon chain (propenyl vs. butenyl)
    total_protons = sum(signal['integration'] for signal in nmr_data.values())
    
    expected_protons_propenyl = 8
    expected_protons_butenyl = 10
    
    derived_chain = None
    if total_protons == expected_protons_propenyl:
        derived_chain = "propenyl"
    elif total_protons == expected_protons_butenyl:
        derived_chain = "butenyl"
    else:
        return f"Incorrect. The total proton count from the NMR data is {total_protons}, which does not match propenyl (8H) or butenyl (10H) acetate."

    # Step 2: Check the J-coupling constant to determine stereochemistry (cis vs. trans)
    # The vinylic coupling constant is given for the signal at 7.0 ppm.
    j_coupling = nmr_data['signal_1']['J_Hz']
    
    # Typical ranges for vinylic coupling
    cis_j_range = (6, 12)
    trans_j_range = (12, 18)
    
    derived_stereochem = None
    if trans_j_range[0] <= j_coupling <= trans_j_range[1]:
        derived_stereochem = "Trans"
    elif cis_j_range[0] <= j_coupling <= cis_j_range[1]:
        derived_stereochem = "Cis"
    else:
        return f"Incorrect. The J-coupling constant is {j_coupling} Hz, which is outside the typical ranges for both cis ({cis_j_range} Hz) and trans ({trans_j_range} Hz) configurations."

    # Step 3: Combine the findings to identify the correct compound
    correct_compound_name = f"{derived_stereochem}-{derived_chain} acetate"

    # Step 4: Validate the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<X>>>' but got '{llm_answer_text}'."
        
    llm_choice_letter = match.group(1)
    llm_choice_name = options.get(llm_choice_letter)

    if not llm_choice_name:
        return f"Invalid option '{llm_choice_letter}' provided in the answer."

    if llm_choice_name.lower() == correct_compound_name.lower():
        return "Correct"
    else:
        reason = []
        # Check chain type
        if derived_chain not in llm_choice_name.lower():
            reason.append(f"the total proton count is {total_protons}, which indicates a '{derived_chain}' structure, not '{'butenyl' if 'butenyl' in llm_choice_name.lower() else 'propenyl'}'")
        # Check stereochemistry
        if derived_stereochem.lower() not in llm_choice_name.lower():
            reason.append(f"the J-coupling constant is {j_coupling} Hz, which indicates a '{derived_stereochem}' configuration, not '{'Cis' if 'cis' in llm_choice_name.lower() else 'Trans'}'")
        
        return (f"Incorrect. The provided answer is {llm_choice_name} ({llm_choice_letter}), but the data points to {correct_compound_name}. "
                f"The answer is wrong because " + " and ".join(reason) + ".")

# Run the check
result = check_answer()
print(result)