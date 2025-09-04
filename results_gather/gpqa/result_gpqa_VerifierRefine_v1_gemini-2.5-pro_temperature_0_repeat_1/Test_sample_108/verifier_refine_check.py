import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the physical constraints for each transition.
    """
    # Define mappings from spectroscopic notation to integer values
    L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

    # The options given in the question
    options = {
        "A": "3D3 + f",
        "B": "3P0 + s",
        "C": "7D1 + p",
        "D": "3S1 + p"
    }
    
    # The answer provided by the LLM
    llm_answer = "B"
    
    # Store the analysis results for each option
    results = {}
    
    for key, transition_str in options.items():
        try:
            # --- 1. Parse the transition string ---
            parts = transition_str.split('+')
            nn_term_str = parts[0].strip()
            x_term_str = parts[1].strip()
            
            # Parse NN final state term symbol: (2S+1)L(J)
            match = re.match(r"(\d+)([A-Z])(\d+)", nn_term_str)
            if not match:
                return f"Error: Could not parse the NN term symbol for option {key}: '{nn_term_str}'"
            
            multiplicity = int(match.group(1))
            L_char = match.group(2)
            J_f = int(match.group(3))
            
            # Calculate S_f from multiplicity 2S+1
            S_f = (multiplicity - 1) / 2
            if S_f != int(S_f):
                 return f"Error: Invalid multiplicity {multiplicity} for option {key}. S_f must be an integer or half-integer."
            S_f = int(S_f)

            # Get L_f from the character
            if L_char not in L_map:
                return f"Error: Unknown L symbol '{L_char}' for option {key}."
            L_f = L_map[L_char]
            
            # Parse particle X orbital angular momentum l_X
            if x_term_str not in l_map:
                return f"Error: Unknown l_X symbol '{x_term_str}' for option {key}."
            l_X = l_map[x_term_str]
            
            # --- 2. Check the three conditions ---
            # Condition 1: Conservation of Angular Momentum (J_f = l_X)
            cond1_ok = (J_f == l_X)
            
            # Condition 2: Conservation of Parity (L_f + l_X is odd)
            cond2_ok = ((L_f + l_X) % 2 == 1)
            
            # Condition 3: Pauli Statistics for T=0 (S_f + L_f is odd)
            cond3_ok = ((S_f + L_f) % 2 == 1)
            
            is_permitted = cond1_ok and cond2_ok and cond3_ok
            
            # Store results and reasons for failure
            results[key] = {
                "permitted": is_permitted,
                "reasons": []
            }
            if not cond1_ok:
                results[key]["reasons"].append(f"Angular momentum not conserved: J_f(NN)={J_f} but l_X={l_X}.")
            if not cond2_ok:
                results[key]["reasons"].append(f"Parity not conserved: L_f + l_X = {L_f} + {l_X} = {L_f + l_X}, which is not odd.")
            if not cond3_ok:
                results[key]["reasons"].append(f"Pauli statistics for T=0 violated: S_f + L_f = {S_f} + {L_f} = {S_f + L_f}, which is not odd.")

        except Exception as e:
            return f"An error occurred while processing option {key}: {e}"

    # --- 3. Final Verification ---
    not_permitted_options = [key for key, res in results.items() if not res["permitted"]]
    
    if len(not_permitted_options) == 0:
        return "The answer is incorrect. The analysis shows all options are permitted, but the question implies one is not."
        
    if len(not_permitted_options) > 1:
        return f"The answer is incorrect. The analysis shows that multiple options are not permitted: {', '.join(not_permitted_options)}."
        
    # Exactly one option is not permitted
    not_permitted_option = not_permitted_options[0]
    
    if not_permitted_option == llm_answer:
        return "Correct"
    else:
        correct_reasoning = results[not_permitted_option]
        return (f"The answer is incorrect. The LLM chose {llm_answer}, but the only non-permitted transition is {not_permitted_option}.\n"
                f"Reason for {not_permitted_option} being not permitted: {' '.join(correct_reasoning['reasons'])}")

# Run the check
result = check_correctness()
print(result)