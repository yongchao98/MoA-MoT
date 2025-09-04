import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the nuclear physics problem.
    The code verifies the provided answer by applying the physical constraints
    to all possible transitions.
    """

    # --- 1. Define helper functions and data ---

    L_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}
    l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

    def parse_term_symbol(symbol):
        """Parses a term symbol like '3P0' into (S, L, J)."""
        match = re.match(r"(\d+)([SPDFG])(\d+)", symbol)
        if not match:
            return None, None, None
        
        multiplicity_str, L_char, J_str = match.groups()
        
        multiplicity = int(multiplicity_str)
        J = int(J_str)
        
        # Calculate S from multiplicity 2S+1
        # This handles integer and half-integer spin possibilities
        S = (multiplicity - 1) / 2
            
        L = L_map.get(L_char)
        
        return S, L, J

    # --- 2. Define the problem parameters and options ---
    
    options = {
        "A": {"nn_final": "3S1", "x_particle": "p"},
        "B": {"nn_final": "3P0", "x_particle": "s"},
        "C": {"nn_final": "7D1", "x_particle": "p"},
        "D": {"nn_final": "3D3", "x_particle": "f"},
    }
    
    llm_answer = "B"

    # --- 3. Apply the physics rules to each option ---

    results = {}
    for option_key, data in options.items():
        nn_symbol = data["nn_final"]
        x_symbol = data["x_particle"]
        
        S_f, L_f, J_f = parse_term_symbol(nn_symbol)
        l_X = l_map.get(x_symbol)
        
        violations = []
        
        # Rule 1: Physicality of NN state (S must be 0 or 1)
        if S_f not in [0, 1]:
            violations.append(f"Physicality rule violated: Final NN spin S_f={S_f} is not possible for a two-nucleon system (must be 0 or 1).")
            
        # Rule 2: Pauli Statistics (S_f + L_f must be odd for T_f=0)
        if (S_f + L_f) % 2 == 0:
            violations.append(f"Pauli statistics rule violated: S_f + L_f = {S_f} + {L_f} = {S_f + L_f}, which is even (must be odd for T_f=0).")
            
        # Rule 3: Parity Conservation (L_f + l_X must be odd)
        if (L_f + l_X) % 2 == 0:
            violations.append(f"Parity conservation rule violated: L_f + l_X = {L_f} + {l_X} = {L_f + l_X}, which is even (must be odd).")
            
        # Rule 4: Angular Momentum Conservation (J_f = l_X)
        if J_f != l_X:
            violations.append(f"Angular momentum conservation rule violated: J_f = {J_f} but l_X = {l_X}.")
            
        results[option_key] = violations

    # --- 4. Evaluate the LLM's answer and reasoning ---

    # Check if the LLM's chosen answer is indeed forbidden.
    if not results[llm_answer]:
        return f"Incorrect. The LLM answered '{llm_answer}', but the code finds this transition to be permitted. Full analysis: {results}"

    # Check if the LLM's reasoning for its answer matches the findings.
    # The LLM's text states B is forbidden due to the Pauli rule.
    llm_reasoning_violation = "Pauli statistics rule"
    found_reason = any(llm_reasoning_violation in v for v in results[llm_answer])
    
    if not found_reason:
        return f"Incorrect. The LLM answered '{llm_answer}' which is a forbidden transition, but the primary reason given in its analysis (violation of Pauli rule) was not found by the code. Found violations: {results[llm_answer]}."

    # Check the LLM's analysis of the ambiguity with option C.
    # The LLM correctly identifies that C is also forbidden due to physicality.
    if "C" in results:
        if not results["C"]:
            return f"Incorrect. The LLM's analysis mentions that option C is also forbidden, but the code finds it to be permitted. Analysis of C: {results['C']}"
        
        llm_reasoning_c = "Physicality rule"
        found_reason_c = any(llm_reasoning_c in v for v in results["C"])
        if not found_reason_c:
            return f"Incorrect. The LLM's analysis of the ambiguity is flawed. It claims C is forbidden due to physicality, but the code finds other/no reasons. Violations for C: {results['C']}"

    # If all checks pass, the LLM's answer and its detailed reasoning are correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)