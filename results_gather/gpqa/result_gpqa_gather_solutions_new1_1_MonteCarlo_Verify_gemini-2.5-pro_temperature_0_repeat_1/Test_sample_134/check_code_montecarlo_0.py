import re

def check_particle_decay_answer():
    """
    Checks the correctness of the answer for the boson decay question.

    The logic is as follows:
    1. Define the mass of the boson X and the resulting mass threshold for the decay products.
    2. Store the standard masses of all fundamental fermions (quarks and leptons).
    3. Determine the set of all kinematically allowed decays by checking which fermions
       have a mass below the threshold.
    4. Parse the proposed answer (Option C) into a set of decays.
    5. Compare the calculated correct set with the answer's set.
    6. Return "Correct" if they match, otherwise return a detailed reason for the error.
    """
    # 1. Define kinematic constraints
    m_X = 6.0  # Mass of boson X in GeV
    fermion_mass_threshold = m_X / 2.0  # m_f must be <= 3.0 GeV

    # 2. Define fermion masses in GeV (approximate values from PDG)
    fermion_data = {
        # symbol: (mass_in_GeV, decay_string_format)
        'e':   (0.000511, 'e^{+}e^{-}'),
        'mu':  (0.1057,   '\mu^{+}\mu^{-}'),
        'tau': (1.777,    '\tau^{+}\tau^{-}'),
        'u':   (0.0022,   'u\bar{u}'),
        'd':   (0.0047,   'd\bar{d}'),
        's':   (0.095,    's\bar{s}'),
        'c':   (1.27,     'c\bar{c}'),
        'b':   (4.18,     'b\bar{b}'),
        't':   (173.0,    't\bar{t}'),
    }

    # 3. Calculate the theoretically correct set of allowed decays
    correct_allowed_decays = set()
    for symbol, (mass, decay_str) in fermion_data.items():
        if mass <= fermion_mass_threshold:
            correct_allowed_decays.add(decay_str)

    # 4. Parse the proposed answer (Option C) into a set
    # The final answer from the LLM is <<<C>>>
    # Option C from the question is: X->c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}
    answer_string = "c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}"
    # Split by comma and remove any leading/trailing whitespace
    answer_decays_set = {decay.strip() for decay in answer_string.split(',')}

    # 5. Compare the correct set with the answer's set
    if correct_allowed_decays == answer_decays_set:
        return "Correct"
    else:
        # Find discrepancies to provide a detailed error message
        missing_decays = correct_allowed_decays - answer_decays_set
        extra_decays = answer_decays_set - correct_allowed_decays

        error_messages = []
        if extra_decays:
            # Find why these are forbidden
            reasons = []
            for decay in sorted(list(extra_decays)):
                for symbol, (mass, decay_str) in fermion_data.items():
                    if decay_str == decay:
                        reasons.append(f"{decay} (fermion mass {mass} GeV > 3 GeV)")
            error_messages.append(f"The answer incorrectly includes forbidden decays: {', '.join(reasons)}.")
        
        if missing_decays:
            error_messages.append(f"The answer omits the following allowed decays: {sorted(list(missing_decays))}.")

        return f"Incorrect. {' '.join(error_messages)}"

# Execute the check and print the result
result = check_particle_decay_answer()
print(result)