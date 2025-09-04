import re

def parse_term_symbol(symbol):
    """Parses a term symbol (2S+1)L(J) into S, L, J quantum numbers."""
    try:
        multiplicity = int(symbol[0])
        S = (multiplicity - 1) / 2
        
        L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3}
        L_char = symbol[1].upper()
        if L_char not in L_map:
            return None, None, None
        L = L_map[L_char]
        
        J = int(symbol[2:])
        return S, L, J
    except (ValueError, IndexError, KeyError):
        return None, None, None

def parse_lx(char):
    """Parses the orbital angular momentum of particle X."""
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
    return l_map.get(char.lower())

def check_transition(transition_str):
    """
    Checks a transition against the physical rules.
    Returns a list of reasons why the transition is not permitted.
    An empty list means the transition is permitted.
    """
    match = re.match(r"1S0 -> (\w+) \+ (\w)", transition_str)
    if not match:
        return ["Invalid format"]
        
    final_nn_symbol = match.group(1)
    particle_x_char = match.group(2)

    S_f, L_f, J_f = parse_term_symbol(final_nn_symbol)
    l_x = parse_lx(particle_x_char)

    if S_f is None or l_x is None:
        return [f"Could not parse the transition string: {transition_str}"]

    violations = []

    # Rule 1: Physicality of the Two-Nucleon State
    # For two nucleons (spin-1/2), total spin S can only be 0 or 1.
    if S_f not in [0, 1]:
        violations.append(f"Unphysical State: Final NN state '{final_nn_symbol}' has spin S={S_f}, which is impossible for two nucleons.")

    # Rule 2: Pauli Statistics for the Final State (T=0)
    # S_f + L_f must be odd.
    if (S_f + L_f) % 2 == 0:
        violations.append(f"Pauli Principle Violation: For T=0, S_f+L_f must be odd, but it is {int(S_f + L_f)} (even).")

    # Rule 3: Conservation of Parity
    # L_f + l_x must be odd.
    if (L_f + l_x) % 2 == 0:
        violations.append(f"Parity Conservation Violation: L_f+l_x must be odd, but it is {int(L_f + l_x)} (even).")

    # Rule 4: Conservation of Total Angular Momentum
    # J_f must equal l_x.
    if J_f != l_x:
        violations.append(f"Angular Momentum Conservation Violation: J_f ({J_f}) must equal l_x ({l_x}).")
        
    return violations

def check_correctness():
    """
    Checks the correctness of the LLM's answer.
    """
    llm_answer = "A"
    
    options = {
        "A": "1S0 -> 3P0 + s",
        "B": "1S0 -> 3D3 + f",
        "C": "1S0 -> 3S1 + p",
        "D": "1S0 -> 7D1 + p"
    }

    # Analyze all options
    results = {}
    forbidden_options = {}
    for option_letter, transition in options.items():
        violations = check_transition(transition)
        results[option_letter] = {
            "permitted": not violations,
            "reasons": violations
        }
        if not violations:
            print(f"Option {option_letter} ({transition}) is PERMITTED.")
        else:
            print(f"Option {option_letter} ({transition}) is NOT PERMITTED.")
            for reason in violations:
                print(f"  - {reason}")
            forbidden_options[option_letter] = violations

    # Check the LLM's answer
    if llm_answer not in forbidden_options:
        return f"Incorrect. The provided answer '{llm_answer}' corresponds to a transition that is permitted by all rules."

    # The LLM correctly identified a forbidden transition.
    # The LLM's reasoning for choosing A over D (both are forbidden) is a matter of interpretation,
    # but the choice itself is valid since the question asks for one non-permitted wave.
    return "Correct"

# Run the check
result = check_correctness()
print(f"\nFinal assessment of the LLM's answer: {result}")
