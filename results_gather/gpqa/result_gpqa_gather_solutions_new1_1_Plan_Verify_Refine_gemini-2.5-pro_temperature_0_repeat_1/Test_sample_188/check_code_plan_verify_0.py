import json

def check_answer():
    """
    Checks the correctness of the provided answer for the physics question.
    """
    # Define the question's options and the provided answer
    question_options = {
        "A": "Magnon",
        "B": "Skyrmion",
        "C": "Pion",
        "D": "Phonon"
    }
    llm_answer_key = "B"

    # Knowledge base based on established physics principles (Goldstone's Theorem)
    # A particle is "associated with SSB" in the context of this question if it is a Goldstone boson.
    knowledge_base = {
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson resulting from the spontaneous breaking of rotational symmetry in ferromagnets."
        },
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "A topological soliton, a stable, particle-like field configuration whose stability comes from topology, not from being a Goldstone boson."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "description": "A pseudo-Goldstone boson from the spontaneous breaking of approximate chiral symmetry in QCD."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson resulting from the spontaneous breaking of continuous translational symmetry in a crystal lattice."
        }
    }

    # Determine the correct answer from the knowledge base
    correct_key = None
    for key, particle_name in question_options.items():
        if not knowledge_base[particle_name]["is_goldstone_boson"]:
            correct_key = key
            break

    # Check if the LLM's answer matches the correct answer
    if llm_answer_key == correct_key:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrectness
        llm_particle_name = question_options.get(llm_answer_key, "an invalid option")
        correct_particle_name = question_options.get(correct_key, "an unknown particle")

        if llm_particle_name == "an invalid option":
            return f"Incorrect. The provided answer key '{llm_answer_key}' is not a valid option."

        llm_particle_info = knowledge_base[llm_particle_name]
        correct_particle_info = knowledge_base[correct_particle_name]

        reason = (
            f"Incorrect. The provided answer is {llm_answer_key} ({llm_particle_name}). "
            f"The question asks which particle is NOT a Goldstone boson resulting from spontaneously-broken symmetry. "
            f"However, {llm_particle_name} is a Goldstone (or pseudo-Goldstone) boson. "
            f"Specifically, it is: {llm_particle_info['description']}. "
            f"The correct answer is {correct_key} ({correct_particle_name}), which is a topological soliton, not a Goldstone boson."
        )
        return reason

# Run the check and print the result
result = check_answer()
print(result)