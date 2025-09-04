import collections

def check_particle_decay_answer():
    """
    Checks the correctness of the selected answer for the particle decay question.

    The core physics principle is that for a decay X -> f(bar)f to be kinematically
    allowed, the mass of the initial particle (m_X) must be greater than or equal
    to the sum of the masses of the final particles (2 * m_f).
    
    Given m_X = 6 GeV, the condition is m_f <= 3 GeV.
    """
    
    # 1. Define the physical constraints and the given answer
    m_X = 6.0  # Mass of boson X in GeV
    final_answer_key = 'C'

    # 2. Define the masses of the fundamental fermions in GeV
    # Using standard approximate values consistent with the provided analysis.
    fermion_masses = {
        'e': 0.000511,    # electron
        'mu': 0.1057,     # muon
        'tau': 1.777,     # tau lepton
        'u': 0.0022,      # up quark
        'd': 0.0047,      # down quark
        's': 0.095,       # strange quark
        'c': 1.27,        # charm quark
        'b': 4.18,        # bottom quark
        't': 173.0,       # top quark
    }

    # 3. Determine the theoretically correct set of allowed decays
    # A decay X -> f(bar)f is allowed if m_f <= m_X / 2
    correctly_allowed_decays = set()
    for fermion_symbol, mass in fermion_masses.items():
        if mass <= m_X / 2:
            correctly_allowed_decays.add(fermion_symbol)

    # 4. Define the decay products listed in each multiple-choice option
    # These are parsed from the question text.
    # Note: 'mu' for muon, 'tau' for tau lepton, 'e' for electron.
    options = {
        'A': {'c', 's', 'u', 'd', 't', 'tau', 'mu', 'e'},
        'B': {'b', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'C': {'c', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'D': {'b', 's', 'u', 'd', 'tau', 'e'}
    }

    # 5. Get the set of decays from the provided final answer
    answer_decays = options.get(final_answer_key)
    if answer_decays is None:
        return f"Error: The final answer key '{final_answer_key}' is not a valid option."

    # 6. Compare the theoretically correct set with the answer's set
    if correctly_allowed_decays == answer_decays:
        return "Correct"
    else:
        # 7. If they don't match, generate a detailed reason for the error
        
        # Decays included in the answer but are kinematically forbidden
        forbidden_included = answer_decays - correctly_allowed_decays
        
        # Decays that are kinematically allowed but were omitted from the answer
        allowed_omitted = correctly_allowed_decays - answer_decays
        
        error_messages = []
        if forbidden_included:
            # Sort for consistent output
            sorted_forbidden = sorted(list(forbidden_included))
            # Provide details on why they are forbidden
            details = [f"{f} (m_f={fermion_masses[f]} GeV)" for f in sorted_forbidden]
            error_messages.append(
                f"The answer incorrectly includes forbidden decays: {details}. "
                f"A decay is forbidden if the fermion's mass is > {m_X/2} GeV."
            )
        
        if allowed_omitted:
            # Sort for consistent output
            sorted_omitted = sorted(list(allowed_omitted))
            # Provide details on why they should be included
            details = [f"{f} (m_f={fermion_masses[f]} GeV)" for f in sorted_omitted]
            error_messages.append(
                f"The answer omits kinematically allowed decays: {details}. "
                f"A decay is allowed if the fermion's mass is <= {m_X/2} GeV."
            )
            
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_particle_decay_answer()
print(result)