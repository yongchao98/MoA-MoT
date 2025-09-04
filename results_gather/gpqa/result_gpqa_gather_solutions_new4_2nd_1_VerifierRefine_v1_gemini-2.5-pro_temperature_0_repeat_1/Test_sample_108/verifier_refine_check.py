import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by applying the selection rules
    of a nuclear transition.
    """

    def parse_term_symbol(symbol):
        """Parses a term symbol like '3D3' into S, L, J."""
        match = re.match(r'(\d+)([A-Z])(\d+)', symbol)
        if not match:
            raise ValueError(f"Invalid term symbol format: {symbol}")
        
        multiplicity_str, L_char, J_str = match.groups()
        
        multiplicity = int(multiplicity_str)
        S = (multiplicity - 1) / 2
        
        L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
        L = L_map[L_char]
        
        J = int(J_str)
        
        return S, L, J

    def get_l_from_wave(wave_char):
        """Converts a wave character like 'p' to its l value."""
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
        return l_map[wave_char]

    def check_transition(final_nn_state, particle_wave):
        """
        Checks a transition against the four rules.
        Returns a list of reasons for failure, or an empty list if permitted.
        """
        failures = []
        
        try:
            S_f, L_f, J_f = parse_term_symbol(final_nn_state)
            l_X = get_l_from_wave(particle_wave)
        except (ValueError, KeyError) as e:
            return [f"Invalid input: {e}"]

        # Rule 1: Physicality of the NN State (S_f must be 0 or 1)
        if S_f not in [0, 1]:
            failures.append(f"Physicality rule violated: Final NN spin S_f={S_f} is not 0 or 1.")

        # Rule 2: Conservation of Total Angular Momentum (J_f = l_X)
        if J_f != l_X:
            failures.append(f"Angular momentum conservation violated: J_f={J_f} but l_X={l_X}.")

        # Rule 3: Conservation of Parity (L_f + l_X must be odd)
        if (L_f + l_X) % 2 == 0:
            failures.append(f"Parity conservation violated: L_f + l_X = {L_f} + {l_X} = {L_f + l_X}, which is not odd.")

        # Rule 4: Pauli Statistics for the Final State (S_f + L_f must be odd for T=0)
        if (S_f + L_f) % 2 == 0:
            failures.append(f"Pauli statistics rule violated: S_f + L_f = {S_f} + {L_f} = {S_f + L_f}, which is not odd for a T=0 final state.")
            
        return failures

    # The options as presented in the question
    options = {
        "A": ("3D3", "f"),
        "B": ("3S1", "p"),
        "C": ("7D1", "p"),
        "D": ("3P0", "s"),
    }
    
    llm_answer = "D"
    
    forbidden_options = {}
    permitted_options = []

    for label, (nn_state, particle_wave) in options.items():
        failures = check_transition(nn_state, particle_wave)
        if not failures:
            permitted_options.append(label)
        else:
            forbidden_options[label] = failures

    # --- Verification Logic ---
    
    # 1. Check if the LLM's answer is indeed forbidden.
    if llm_answer not in forbidden_options:
        return f"Incorrect. The provided answer '{llm_answer}' is a permitted transition. The code found {permitted_options} to be permitted."

    # 2. Check the reasoning for the LLM's answer (D).
    # The analysis states D is forbidden by the Pauli rule ONLY.
    failures_for_D = forbidden_options.get(llm_answer, [])
    pauli_failure_for_D = any("Pauli statistics" in f for f in failures_for_D)
    if not (pauli_failure_for_D and len(failures_for_D) == 1):
        return f"Incorrect reasoning for answer '{llm_answer}'. The code found the following failures: {failures_for_D}. The expected failure was a single violation of the Pauli statistics rule."

    # 3. Check the reasoning for the other forbidden option (C).
    # The analysis states C is forbidden by the physicality rule ONLY.
    if "C" not in forbidden_options:
        return f"Incorrect. The analysis missed that option 'C' is also forbidden. Forbidden options found by code: {list(forbidden_options.keys())}."
    
    failures_for_C = forbidden_options.get("C", [])
    physicality_failure_for_C = any("Physicality rule" in f for f in failures_for_C)
    # For C, the other rules (J, Parity, Pauli) are satisfied.
    if not (physicality_failure_for_C and len(failures_for_C) == 1):
         return f"Incorrect reasoning for option 'C'. The code found the following failures: {failures_for_C}. The expected failure was a single violation of the physicality rule."

    # 4. Check that the permitted options are correct.
    if set(permitted_options) != {"A", "B"}:
        return f"Incorrect. The analysis of permitted transitions is wrong. Code found {permitted_options} to be permitted, but expected ['A', 'B']."

    # If all checks pass, the LLM's answer and detailed reasoning are correct.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)