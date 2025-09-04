import sys
from io import StringIO

def check_answer():
    """
    This function checks the correctness of the provided answer for the 1H NMR problem.
    It verifies the reasoning by checking proton count and stereochemistry based on J-coupling.
    """
    
    # --- Data Definition ---

    # 1. Experimental Data from the question
    experimental_data = {
        'signals': [
            {'ppm': 7.0, 'integration': 1, 'multiplicity': 'd', 'J_Hz': 16.0},
            {'ppm': 5.5, 'integration': 1, 'multiplicity': 'dq'},
            {'ppm': 2.1, 'integration': 3, 'multiplicity': 's'},
            {'ppm': 1.6, 'integration': 3, 'multiplicity': 'd'}
        ]
    }
    
    # 2. Properties of the candidate molecules
    candidates = {
        'A': {'name': 'Trans-propenyl acetate', 'total_protons': 8, 'stereochemistry': 'trans'},
        'B': {'name': 'Cis-butenyl acetate', 'total_protons': 10, 'stereochemistry': 'cis'},
        'C': {'name': 'Cis-propenyl acetate', 'total_protons': 8, 'stereochemistry': 'cis'},
        'D': {'name': 'Trans-butenyl acetate', 'total_protons': 10, 'stereochemistry': 'trans'}
    }
    
    # 3. The answer provided by the LLM
    llm_answer_key = 'A'
    
    # --- Verification Logic ---
    
    # Step 1: Check total proton count
    experimental_protons = sum(s['integration'] for s in experimental_data['signals'])
    
    if experimental_protons != 8:
        return f"Incorrect Pre-analysis: The total integration from the experimental data is {experimental_protons}, not 8. Please check the signal integrations."

    # Get the properties of the selected answer
    selected_candidate = candidates.get(llm_answer_key)
    if not selected_candidate:
        return f"Invalid Answer Key: The key '{llm_answer_key}' does not correspond to any of the choices."

    # Check if the selected answer's proton count matches the experimental data
    if selected_candidate['total_protons'] != experimental_protons:
        return (f"Incorrect based on proton count. The experimental data shows {experimental_protons} protons. "
                f"The selected answer, {selected_candidate['name']}, has {selected_candidate['total_protons']} protons.")

    # Verify that this check correctly eliminates other candidates
    eliminated_by_protons = [key for key, props in candidates.items() if props['total_protons'] != experimental_protons]
    if sorted(eliminated_by_protons) != ['B', 'D']:
        return (f"Reasoning Error: The proton count check should eliminate candidates B and D (10 protons), "
                f"but the code found that {eliminated_by_protons} were eliminated.")

    # Step 2: Check stereochemistry via J-coupling
    # Typical J-coupling for vinylic protons: trans = 12-18 Hz, cis = 6-12 Hz
    vinylic_coupling_signal = next((s for s in experimental_data['signals'] if s.get('J_Hz')), None)
    
    if not vinylic_coupling_signal:
        return "Constraint Check Failed: Could not find a signal with a J-coupling constant in the experimental data to determine stereochemistry."

    j_value = vinylic_coupling_signal['J_Hz']
    
    determined_stereochemistry = ''
    if 12 <= j_value <= 18:
        determined_stereochemistry = 'trans'
    elif 6 <= j_value < 12:
        determined_stereochemistry = 'cis'
    else:
        return f"Constraint Check Failed: The J-coupling value of {j_value} Hz is unusual and does not clearly fit typical cis or trans ranges."

    # Check if the selected answer's stereochemistry matches the data
    if selected_candidate['stereochemistry'] != determined_stereochemistry:
        return (f"Incorrect based on stereochemistry. The J-coupling of {j_value} Hz indicates a '{determined_stereochemistry}' configuration, "
                f"but the selected answer is '{selected_candidate['name']}' which is '{selected_candidate['stereochemistry']}'.")

    # Verify that this check correctly eliminates the remaining incorrect candidate
    remaining_candidates = {k: v for k, v in candidates.items() if k not in eliminated_by_protons}
    eliminated_by_j_coupling = [key for key, props in remaining_candidates.items() if props['stereochemistry'] != determined_stereochemistry]
    
    if 'C' not in eliminated_by_j_coupling:
         return (f"Reasoning Error: The J-coupling check should eliminate candidate C (cis), but it was not eliminated.")

    # Step 3: Final Verification
    # At this point, only the correct answer should remain.
    final_candidate_key = [k for k, v in candidates.items() if k not in eliminated_by_protons and k not in eliminated_by_j_coupling]
    
    if len(final_candidate_key) != 1 or final_candidate_key[0] != llm_answer_key:
        return (f"Reasoning Error: After applying all constraints, the remaining candidate should be '{llm_answer_key}', "
                f"but the analysis resulted in: {final_candidate_key}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)