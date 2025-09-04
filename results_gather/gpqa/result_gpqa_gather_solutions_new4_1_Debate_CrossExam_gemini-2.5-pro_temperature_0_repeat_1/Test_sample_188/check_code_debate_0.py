import collections

def check_answer():
    """
    This function checks the correctness of the LLM's answer to a physics question.
    It uses a predefined knowledge base about effective particles and their relation to
    spontaneous symmetry breaking (SSB).
    """

    # Knowledge base: Defines the relationship of each particle to SSB.
    # The key distinction is whether the particle is a Goldstone boson (or pseudo-GB),
    # which is a direct consequence of a spontaneously broken continuous symmetry.
    particle_definitions = {
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "a topological soliton, whose stability comes from topology, not Goldstone's theorem."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "the Goldstone boson of spontaneously broken translational symmetry in a crystal."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "description": "the pseudo-Goldstone boson of spontaneously broken approximate chiral symmetry in QCD."
        },
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "the Goldstone boson of spontaneously broken rotational symmetry in a ferromagnet."
        }
    }

    # The options as presented in the question.
    options = {
        "A": "Skyrmion",
        "B": "Phonon",
        "C": "Pion",
        "D": "Magnon"
    }

    # The final answer provided by the LLM analysis.
    llm_answer_key = "A"

    # Determine the correct answer based on the knowledge base.
    # The question asks which particle is NOT a Goldstone boson.
    correct_particle_name = None
    for particle, properties in particle_definitions.items():
        if not properties["is_goldstone_boson"]:
            correct_particle_name = particle
            break
    
    # Find the key (A, B, C, or D) corresponding to the correct particle name.
    correct_key = None
    for key, value in options.items():
        if value == correct_particle_name:
            correct_key = key
            break

    # Check if the LLM's answer matches the correct answer.
    if llm_answer_key == correct_key:
        return "Correct"
    else:
        llm_particle_name = options.get(llm_answer_key, "Invalid Option")
        correct_particle_description = particle_definitions[correct_particle_name]["description"]
        llm_particle_description = particle_definitions[llm_particle_name]["description"]

        reason = (
            f"Incorrect. The provided answer is {llm_answer_key} ({llm_particle_name}), but the correct answer is {correct_key} ({correct_particle_name}).\n"
            f"The question asks to identify the particle that is not a Goldstone boson resulting from spontaneous symmetry breaking.\n"
            f"- The correct answer is the Skyrmion because it is {correct_particle_description}\n"
            f"- The provided answer, {llm_particle_name}, is incorrect because it is {llm_particle_description}\n"
            f"- Phonons, Pions, and Magnons are all (pseudo-)Goldstone bosons, which are the direct particle manifestations of spontaneously broken continuous symmetries."
        )
        return reason

# Execute the check and print the result.
result = check_answer()
print(result)