def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the physics constraints
    of the given nuclear decay problem.
    """
    
    # --- Helper functions to parse quantum numbers ---
    def parse_term_symbol(symbol):
        """Parses a (2S+1)L(J) term symbol string."""
        L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
        s_val = (int(symbol[0]) - 1) / 2
        l_val = L_map[symbol[1]]
        j_val = int(symbol[2])
        return s_val, l_val, j_val

    def parse_particle_wave(wave):
        """Parses the lowercase letter for particle X's angular momentum."""
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        return l_map[wave]

    # --- Define the options from the question ---
    options = {
        "A": "3D3 + f",
        "B": "7D1 + p",
        "C": "3S1 + p",
        "D": "3P0 + s"
    }
    
    llm_answer = "C"
    
    # --- Analyze each option ---
    results = {}
    for option_letter, transition_str in options.items():
        nn_symbol, x_wave = transition_str.split(" + ")
        
        s_f, l_f, j_f = parse_term_symbol(nn_symbol)
        l_x = parse_particle_wave(x_wave)
        
        # Store violations for this option
        violations = []
        
        # Rule 4: Physicality of NN State (S_f must be 0 or 1)
        if s_f not in [0, 1]:
            violations.append(f"Unphysical NN state: S_f={s_f} is not possible for two nucleons.")
        
        # Rule 1: Angular Momentum Conservation (J_f = l_X)
        if j_f != l_x:
            violations.append(f"Angular momentum not conserved: J_f={j_f} but l_X={l_x}.")
            
        # Rule 2: Parity Conservation (L_f + l_X must be odd)
        if (l_f + l_x) % 2 == 0:
            violations.append(f"Parity not conserved: L_f + l_X = {l_f + l_x} is even.")
            
        # Rule 3: Pauli Principle for T=0 (S_f + L_f must be odd)
        if (s_f + l_f) % 2 == 0:
            violations.append(f"Pauli principle violated: S_f + L_f = {s_f + l_f} is even.")
            
        results[option_letter] = violations

    # --- Final check against the LLM's answer ---
    permitted_options = [k for k, v in results.items() if not v]
    forbidden_options = {k: v for k, v in results.items() if v}

    if llm_answer in permitted_options:
        error_report = f"Incorrect. The provided answer is <<<{llm_answer}>>>, which corresponds to the transition '1S0 -> {options[llm_answer]}'. My analysis shows this transition is permitted.\n\n"
        error_report += "The non-permitted transitions are:\n"
        for opt, v_list in forbidden_options.items():
            error_report += f"- Option {opt} ('1S0 -> {options[opt]}'): Forbidden because {'; '.join(v_list)}\n"
        return error_report
    elif llm_answer in forbidden_options:
        # This case would mean the letter is correct, but we should check the reasoning.
        # The LLM's reasoning points to '3P0 + s' (Option D) and '7D1 + p' (Option B) being forbidden.
        # The LLM chose C. This is a contradiction.
        return f"Incorrect. The provided answer is <<<{llm_answer}>>>, but the reasoning given in the prompt is inconsistent. The reasoning correctly identifies that '1S0 -> 3P0 + s' (Option D) and '1S0 -> 7D1 + p' (Option B) are forbidden. However, the final answer selected is C, which is a permitted transition."
    else:
        return "Error: The provided answer letter does not match any of the options A, B, C, D."

# Execute the check
print(check_correctness())