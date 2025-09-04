import re

def parse_term_symbol(symbol):
    """
    Parses a term symbol like '3P0' into its S, L, and J quantum numbers.
    Returns (S, L, J).
    """
    try:
        match = re.match(r'(\d+)([SPDFG])(\d+)', symbol)
        if not match:
            return None, None, None
            
        multiplicity = int(match.group(1))
        l_char = match.group(2)
        j_val = int(match.group(3))
        
        s_val = (multiplicity - 1) / 2
        
        l_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
        l_val = l_map[l_char]
        
        return s_val, l_val, j_val
    except (KeyError, ValueError, IndexError):
        return None, None, None

def parse_particle_wave(wave):
    """
    Parses a particle wave character like 'p' into its orbital angular momentum.
    Returns l_x.
    """
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    return l_map.get(wave, None)

def check_transition(final_nn_state, particle_wave):
    """
    Checks a given transition against the four rules.
    Returns a list of reasons why it's forbidden, or an empty list if permitted.
    """
    reasons = []
    
    S_f, L_f, J_f = parse_term_symbol(final_nn_state)
    l_x = parse_particle_wave(particle_wave)
    
    if S_f is None or l_x is None:
        return [f"Invalid state notation: {final_nn_state} or {particle_wave}"]

    # Rule 4: Physicality of the NN state (two nucleons are spin-1/2)
    if S_f not in [0, 1]:
        reasons.append(f"Unphysical NN state: Total spin S_f={S_f} is not possible for two nucleons (must be 0 or 1).")
        # Even if unphysical, we can check the other rules for completeness.

    # Rule 1: Conservation of Total Angular Momentum (J_f = l_x)
    if J_f != l_x:
        reasons.append(f"Violates J conservation: J_f={J_f} but l_x={l_x}.")

    # Rule 2: Conservation of Parity (L_f + l_x must be odd)
    if (L_f + l_x) % 2 == 0:
        reasons.append(f"Violates parity conservation: L_f + l_x = {L_f + l_x} is even, but must be odd.")

    # Rule 3: Pauli Principle for T=0 (S_f + L_f must be odd)
    if (S_f + L_f) % 2 == 0:
        reasons.append(f"Violates Pauli principle for T=0: S_f + L_f = {S_f + L_f} is even, but must be odd.")
        
    return reasons

def check_correctness():
    """
    Main function to check the LLM's answer and reasoning.
    """
    llm_answer = "B"
    
    options = {
        "A": ("7D1", "p"),
        "B": ("3P0", "s"),
        "C": ("3D3", "f"),
        "D": ("3S1", "p")
    }
    
    results = {}
    for option, (nn_state, particle) in options.items():
        reasons = check_transition(nn_state, particle)
        results[option] = {
            "permitted": len(reasons) == 0,
            "reasons": reasons
        }
        
    # 1. Check if the LLM's chosen answer is indeed not permitted.
    if results[llm_answer]["permitted"]:
        return f"Incorrect. The provided answer '{llm_answer}' is a permitted transition according to the rules. My analysis: {results}"
        
    # 2. Check if the LLM's reasoning for its answer is correct.
    # The LLM states B is forbidden due to violating the Pauli principle.
    llm_reason_for_b = "Violates the Pauli principle"
    found_reasons_for_b = results[llm_answer]["reasons"]
    
    pauli_violation_found = any("Pauli principle" in r for r in found_reasons_for_b)
    
    if not pauli_violation_found:
        return f"Incorrect. The provided answer '{llm_answer}' is indeed not permitted, but the stated reason is wrong. The LLM claims it violates the Pauli principle, but my analysis found these reasons: {found_reasons_for_b}."

    # 3. Check the LLM's analysis of the other forbidden state (A).
    # The LLM correctly identifies that A is also forbidden because it's unphysical.
    if results["A"]["permitted"]:
        return f"Incorrect. The LLM correctly identifies that B is forbidden, but misses that A is also forbidden. My analysis shows A is forbidden for reasons: {results['A']['reasons']}"
    
    unphysical_violation_found = any("Unphysical" in r for r in results["A"]["reasons"])
    if not unphysical_violation_found:
        return f"Incorrect. The LLM's reasoning about option A being unphysical is not supported by my analysis. My reasons for A: {results['A']['reasons']}."

    # 4. Check the LLM's analysis of the permitted states (C and D).
    if not results["C"]["permitted"]:
        return f"Incorrect. The LLM's analysis is flawed. My analysis shows that C is not permitted for reasons: {results['C']['reasons']}."
    if not results["D"]["permitted"]:
        return f"Incorrect. The LLM's analysis is flawed. My analysis shows that D is not permitted for reasons: {results['D']['reasons']}."

    # 5. If all checks pass, the LLM's answer and its detailed reasoning are correct.
    return "Correct"

# Run the check
print(check_correctness())