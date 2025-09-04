import re

def check_correctness_of_answer():
    """
    Checks the correctness of the final answer for the given physics question.

    The function uses a knowledge base to determine the correct answer based on
    the principles of Spontaneous Symmetry Breaking and Goldstone's Theorem.
    It then compares this correct answer with the provided one.
    """

    # Knowledge base defining the nature of each particle in relation to SSB.
    # The question asks which particle is NOT a Goldstone boson.
    particle_info = {
        'Pion': {
            'is_goldstone_boson': True,
            'description': "A pseudo-Goldstone boson from the spontaneous breaking of approximate chiral symmetry in QCD."
        },
        'Magnon': {
            'is_goldstone_boson': True,
            'description': "A Goldstone boson from the spontaneous breaking of rotational symmetry in a ferromagnet."
        },
        'Phonon': {
            'is_goldstone_boson': True,
            'description': "A Goldstone boson from the spontaneous breaking of continuous translational symmetry in a crystal lattice."
        },
        'Skyrmion': {
            'is_goldstone_boson': False,
            'description': "A topological soliton, whose stability comes from topology, not Goldstone's theorem. It is not a Goldstone boson."
        }
    }

    # The options as presented in the question.
    options = {
        'A': 'Pion',
        'B': 'Magnon',
        'C': 'Phonon',
        'D': 'Skyrmion'
    }

    # The final answer provided by the LLM.
    final_answer_text = "<<<D>>>"

    # Step 1: Determine the correct answer from the knowledge base.
    # The correct particle is the one that is NOT a Goldstone boson.
    correct_particle_name = None
    for particle, info in particle_info.items():
        if not info['is_goldstone_boson']:
            correct_particle_name = particle
            break

    correct_option_letter = None
    for letter, name in options.items():
        if name == correct_particle_name:
            correct_option_letter = letter
            break

    if not correct_option_letter:
        # This would indicate an issue with the checker's internal logic or data.
        return "Error in checker: Could not map the correct particle to an option."

    # Step 2: Extract the given answer from the text.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Incorrect format: The final answer '{final_answer_text}' is not in the format '<<<X>>>'."

    provided_answer_letter = match.group(1)

    # Step 3: Compare the provided answer with the correct answer.
    if provided_answer_letter == correct_option_letter:
        return "Correct"
    else:
        provided_particle_name = options.get(provided_answer_letter, "Unknown")
        reason = (
            f"The provided answer '{provided_answer_letter}' ({provided_particle_name}) is incorrect.\n"
            f"The correct answer is '{correct_option_letter}' ({correct_particle_name}).\n\n"
            f"Reasoning:\n"
            f"The question asks which particle is not a Goldstone boson, which is a particle that arises from spontaneously broken continuous symmetry.\n"
            f"- Pion, Magnon, and Phonon are all examples of Goldstone (or pseudo-Goldstone) bosons.\n"
            f"- A Skyrmion is a topological soliton, not a Goldstone boson. Its stability is due to its topology.\n"
            f"Therefore, Skyrmion is the particle that is not associated with a spontaneously-broken symmetry in this way."
        )
        return reason

# Execute the check and print the result.
result = check_correctness_of_answer()
print(result)