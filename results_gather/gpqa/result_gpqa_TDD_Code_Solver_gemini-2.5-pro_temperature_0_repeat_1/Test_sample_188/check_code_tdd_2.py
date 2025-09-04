def check_answer():
    """
    Checks the correctness of the answer to the particle physics question.
    """
    # Database of particle properties related to Spontaneously Broken Symmetry (SBS)
    # A particle is associated with SBS if it's a Goldstone or pseudo-Goldstone boson.
    particle_info = {
        'Phonon': {
            'is_sbs_associated': True,
            'reason': "It is the Goldstone boson of spontaneously broken translational symmetry in a crystal."
        },
        'Magnon': {
            'is_sbs_associated': True,
            'reason': "It is the Goldstone boson of spontaneously broken rotational symmetry in a ferromagnet."
        },
        'Pion': {
            'is_sbs_associated': True,
            'reason': "It is the pseudo-Goldstone boson of spontaneously broken approximate chiral symmetry in QCD."
        },
        'Skyrmion': {
            'is_sbs_associated': False,
            'reason': "It is a topological soliton, not a Goldstone boson. Its stability is topological."
        }
    }

    # The question's options and the given answer
    options = {
        'A': 'Phonon',
        'B': 'Magnon',
        'C': 'Pion',
        'D': 'Skyrmion'
    }
    given_answer_key = 'D'

    # 1. Identify the particle selected by the given answer.
    selected_particle = options.get(given_answer_key)
    if not selected_particle:
        return f"Invalid answer key '{given_answer_key}'. Options are A, B, C, D."

    # 2. Check if the selected particle satisfies the condition (is NOT associated with SBS).
    if particle_info[selected_particle]['is_sbs_associated']:
        return (f"Incorrect. The selected particle '{selected_particle}' (Answer {given_answer_key}) "
                f"is associated with spontaneously-broken symmetry. "
                f"Reason: {particle_info[selected_particle]['reason']}")

    # 3. Verify that all other options do NOT satisfy the condition (i.e., they ARE associated with SBS).
    # This ensures the answer is uniquely correct among the choices.
    for key, particle_name in options.items():
        if key != given_answer_key:
            if not particle_info[particle_name]['is_sbs_associated']:
                return (f"Incorrect. The answer '{selected_particle}' is not uniquely correct. "
                        f"The particle '{particle_name}' (Option {key}) is also not associated with SBS, "
                        f"which contradicts the premise of the question having a single correct answer.")

    # 4. If the selected particle is the only one not associated with SBS, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)