import re

def check_correctness_of_particle_decay_answer():
    """
    This function checks the correctness of the provided answer for the particle decay problem.

    The logic is as follows:
    1.  Define the physical constants: the mass of the boson X and the masses of all relevant fermions.
    2.  Calculate the kinematic threshold for decay: for a decay X -> f-f_bar, the mass of the fermion f (m_f)
        must be less than or equal to half the mass of the boson X (m_X / 2).
    3.  Determine the theoretically correct set of allowed decays by comparing each fermion's mass to the threshold.
    4.  Parse the multiple-choice options provided in the question to create a mapping from letters (A, B, C, D)
        to the set of decays they represent.
    5.  Extract the letter from the provided final answer (e.g., 'C' from '<<<C>>>').
    6.  Compare the set of decays from the chosen answer option with the theoretically correct set.
    7.  If they match, the answer is correct. Otherwise, provide a detailed reason for the discrepancy,
        specifying which decays were incorrectly included or omitted.
    """
    # 1. Define physical constants
    m_X = 6.0  # Mass of boson X in GeV
    fermion_masses = {
        # Leptons (in GeV)
        'e': 0.000511,
        'mu': 0.1057,
        'tau': 1.777,
        # Quarks (in GeV)
        'u': 0.0022,
        'd': 0.0047,
        's': 0.095,
        'c': 1.27,
        'b': 4.18,
        't': 173.0,
    }

    # 2. Calculate the kinematic threshold
    mass_threshold = m_X / 2.0

    # 3. Determine the theoretically correct set of allowed decays
    correctly_allowed_decays = {
        fermion for fermion, mass in fermion_masses.items() if mass <= mass_threshold
    }

    # 4. Parse the multiple-choice options from the question text
    # Note: The options are taken from the final answer's analysis, which correctly reflects the question.
    options_text = {
        'A': r"b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},e^{+}e^{-}",
        'B': r"b\bar{b},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        'C': r"c\bar{c},s\bar{s},u\bar{u},d\bar{d},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}",
        'D': r"c\bar{c},s\bar{s},u\bar{u},d\bar{d},t\bar{t},\tau^{+}\tau^{-},\mu^{+}\mu^{-},e^{+}e^{-}"
    }

    def parse_decay_string(decay_str):
        """A helper function to parse the decay string into a set of particle names."""
        particles = set()
        # Split by comma and trim whitespace
        parts = [p.strip() for p in decay_str.split(',')]
        for part in parts:
            if part.startswith('b'): particles.add('b')
            elif part.startswith('s'): particles.add('s')
            elif part.startswith('u'): particles.add('u')
            elif part.startswith('d'): particles.add('d')
            elif part.startswith('c'): particles.add('c')
            elif part.startswith('t'): particles.add('t')
            elif part.startswith(r'\tau'): particles.add('tau')
            elif part.startswith(r'\mu'): particles.add('mu')
            elif part.startswith('e'): particles.add('e')
        return particles

    parsed_options = {key: parse_decay_string(value) for key, value in options_text.items()}

    # 5. Extract the letter from the provided final answer
    final_answer_letter = 'C'
    answer_decays = parsed_options.get(final_answer_letter)

    # 6. Compare the answer's decay set with the correct set
    if answer_decays == correctly_allowed_decays:
        return "Correct"
    else:
        # 7. Provide a detailed reason for the error
        error_messages = []
        
        # Check for incorrectly included (forbidden) decays
        extra_decays = answer_decays - correctly_allowed_decays
        if extra_decays:
            for decay in sorted(list(extra_decays)) :
                mass = fermion_masses[decay]
                error_messages.append(
                    f"it incorrectly includes the forbidden decay to '{decay}' (mass ≈ {mass:.2f} GeV > {mass_threshold} GeV)"
                )

        # Check for incorrectly omitted (allowed) decays
        missing_decays = correctly_allowed_decays - answer_decays
        if missing_decays:
            for decay in sorted(list(missing_decays)):
                mass = fermion_masses[decay]
                error_messages.append(
                    f"it incorrectly omits the allowed decay to '{decay}' (mass ≈ {mass:.2f} GeV ≤ {mass_threshold} GeV)"
                )
        
        reason = f"Incorrect. The chosen answer '{final_answer_letter}' is wrong because " + " and ".join(error_messages) + "."
        return reason

# Run the checker
result = check_correctness_of_particle_decay_answer()
print(result)