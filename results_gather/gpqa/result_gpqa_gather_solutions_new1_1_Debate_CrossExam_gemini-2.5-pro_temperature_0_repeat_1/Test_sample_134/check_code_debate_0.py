def check_particle_decay_answer():
    """
    Checks the correctness of the answer for the kinematically allowed decays of boson X.
    """
    # Mass of the boson X in GeV
    m_X = 6.0

    # The kinematic condition is that the fermion mass m_f must be <= m_X / 2
    mass_threshold = m_X / 2.0

    # Approximate masses of fundamental fermions in GeV
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

    # Determine the set of theoretically allowed decays
    allowed_decays = set()
    for fermion, mass in fermion_masses.items():
        if mass <= mass_threshold:
            allowed_decays.add(fermion)

    # The final answer from the LLM to be checked
    llm_answer_choice = 'A'

    # Define the decay products listed in each option as sets of strings
    options = {
        'A': {'c', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'B': {'c', 's', 'u', 'd', 't', 'tau', 'mu', 'e'},
        'C': {'b', 's', 'u', 'd', 'tau', 'mu', 'e'},
        'D': {'b', 's', 'u', 'd', 'tau', 'e'}
    }

    # Get the set of decays corresponding to the LLM's answer
    llm_answer_set = options.get(llm_answer_choice)

    # Compare the theoretically calculated set with the LLM's answer set
    if allowed_decays == llm_answer_set:
        return "Correct"
    else:
        # Find the specific reasons for the mismatch
        missing_decays = allowed_decays - llm_answer_set
        forbidden_included_decays = llm_answer_set - allowed_decays

        reasons = []
        if forbidden_included_decays:
            for f in sorted(list(forbidden_included_decays)):
                reasons.append(
                    f"it incorrectly includes the forbidden decay to '{f}' "
                    f"(since its mass {fermion_masses[f]} GeV > {mass_threshold} GeV)"
                )
        
        if missing_decays:
            for f in sorted(list(missing_decays)):
                reasons.append(
                    f"it omits the allowed decay to '{f}' "
                    f"(since its mass {fermion_masses[f]} GeV <= {mass_threshold} GeV)"
                )
        
        return f"The answer '{llm_answer_choice}' is incorrect because " + ", and ".join(reasons) + "."

# Run the check
result = check_particle_decay_answer()
print(result)