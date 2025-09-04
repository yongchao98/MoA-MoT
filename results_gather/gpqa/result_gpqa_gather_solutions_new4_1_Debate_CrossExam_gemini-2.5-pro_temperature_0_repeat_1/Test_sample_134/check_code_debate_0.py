import re

def check_particle_decay_answer():
    """
    Checks the correctness of the final answer for the particle decay question.

    The core principle is that a boson X with mass m_X can decay into a
    fermion-antifermion pair (f-f_bar) only if the decay is kinematically
    allowed. The condition is m_X >= 2 * m_f, where m_f is the mass of the
    fermion.
    """

    # Mass of the boson X in GeV
    m_X = 6.0

    # Fermion masses in GeV (using values consistent with the provided answers)
    fermion_masses = {
        'e': 0.000511,    # electron
        'mu': 0.1057,     # muon
        'tau': 1.78,      # tau lepton
        'u': 0.0022,      # up quark
        'd': 0.0047,      # down quark
        's': 0.095,       # strange quark
        'c': 1.27,        # charm quark
        'b': 4.18,        # bottom quark
        't': 173.0,       # top quark
    }

    # --- Step 1: Determine the theoretically correct set of allowed decays ---
    
    # The condition is m_f <= m_X / 2
    mass_threshold = m_X / 2
    
    theoretically_allowed_decays = set()
    for fermion, mass in fermion_masses.items():
        if mass <= mass_threshold:
            theoretically_allowed_decays.add(fermion)

    # --- Step 2: Parse the options from the original question ---
    # A) X->c c,s s,u u,d d,tau+ tau-,mu+ mu-,e+ e-
    # B) X->b b,s s,u u,d d,tau+ tau-,e+ e-
    # C) X->b b,s s,u u,d d,tau+ tau-,mu+ mu-,e+ e-
    # D) X->c c,s s,u u,d d,t t,tau+ tau-,mu+ mu-,e+ e-
    
    # Helper to parse decay strings like 'c c' or 'tau+ tau-' into a single identifier
    def parse_decay_string(s):
        s = s.replace('+', '').replace('-', '').replace(' ', '')
        if 'mu' in s: return 'mu'
        if 'tau' in s: return 'tau'
        return s[0]

    options_text = {
        "A": "c c,s s,u u,d d,tau+ tau-,mu+ mu-,e+ e-",
        "B": "b b,s s,u u,d d,tau+ tau-,e+ e-",
        "C": "b b,s s,u u,d d,tau+ tau-,mu+ mu-,e+ e-",
        "D": "c c,s s,u u,d d,t t,tau+ tau-,mu+ mu-,e+ e-"
    }
    
    parsed_options = {}
    for key, text in options_text.items():
        # Split by comma, then parse each part
        decays = [parse_decay_string(part.strip()) for part in text.split(',')]
        parsed_options[key] = set(decays)

    # --- Step 3: Check the provided final answer ---
    final_answer_key = "A"
    llm_answer_set = parsed_options.get(final_answer_key)

    if llm_answer_set is None:
        return f"The final answer key '{final_answer_key}' is not a valid option (A, B, C, or D)."

    # Compare the LLM's answer with the theoretically correct set
    if llm_answer_set == theoretically_allowed_decays:
        return "Correct"
    else:
        # Find the differences to provide a clear reason
        missing_decays = theoretically_allowed_decays - llm_answer_set
        forbidden_included_decays = llm_answer_set - theoretically_allowed_decays
        
        error_messages = []
        if forbidden_included_decays:
            for decay in sorted(list(forbidden_included_decays)):
                mass = fermion_masses[decay]
                error_messages.append(
                    f"The answer incorrectly includes the decay to '{decay}' (mass={mass} GeV). "
                    f"This decay is forbidden because the fermion's mass ({mass} GeV) is greater than the threshold of {mass_threshold} GeV."
                )
        
        if missing_decays:
            for decay in sorted(list(missing_decays)):
                mass = fermion_masses[decay]
                error_messages.append(
                    f"The answer is missing the allowed decay to '{decay}' (mass={mass} GeV). "
                    f"This decay is allowed because the fermion's mass ({mass} GeV) is less than or equal to the threshold of {mass_threshold} GeV."
                )
        
        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result
result = check_particle_decay_answer()
print(result)