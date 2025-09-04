import sys
import io

def check_particle_symmetry_answer():
    """
    This function checks the correctness of the answer to the physics question about
    particles and spontaneously broken symmetry.
    """
    # Define the knowledge base about the particles and their relation to spontaneous symmetry breaking (SSB).
    # The key concept is Goldstone's theorem, which states that SSB of a continuous symmetry
    # gives rise to massless (or low-mass) particles called Goldstone bosons.
    # The question asks which particle is NOT a Goldstone boson.
    knowledge_base = {
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "A magnon is the Goldstone boson from the spontaneous breaking of continuous rotational symmetry in ferromagnets."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "An acoustic phonon is the Goldstone boson from the spontaneous breaking of continuous translational symmetry in a crystal lattice."
        },
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "A skyrmion is a topological soliton, not a Goldstone boson. Its stability is due to its topological properties, not Goldstone's theorem."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "description": "A pion is a pseudo-Goldstone boson from the spontaneous breaking of approximate chiral symmetry in quantum chromodynamics (QCD)."
        }
    }

    # The options provided in the question
    options = {
        "A": "Magnon",
        "B": "Phonon",
        "C": "Skyrmion",
        "D": "Pion"
    }

    # The answer provided by the LLM to be checked
    llm_answer_letter = "C"

    # Determine the correct answer based on the knowledge base
    correct_letter = None
    for letter, particle_name in options.items():
        if not knowledge_base[particle_name]["is_goldstone_boson"]:
            correct_letter = letter
            break  # Assume there is only one correct answer

    # Verify the LLM's answer
    if llm_answer_letter == correct_letter:
        return "Correct"
    else:
        llm_particle_name = options.get(llm_answer_letter)
        if not llm_particle_name:
            return f"Invalid answer: The provided answer '{llm_answer_letter}' is not a valid option."

        correct_particle_name = options[correct_letter]
        
        reason = (
            f"The answer '{llm_answer_letter}' ({llm_particle_name}) is incorrect. "
            f"The question asks which particle is NOT a Goldstone boson resulting from spontaneous symmetry breaking. "
            f"The particle {llm_particle_name} IS considered a Goldstone (or pseudo-Goldstone) boson. "
            f"Details: {knowledge_base[llm_particle_name]['description']}. "
            f"The correct answer is '{correct_letter}' ({correct_particle_name}). "
            f"Details: {knowledge_base[correct_particle_name]['description']}."
        )
        return reason

# Execute the check and print the result
result = check_particle_symmetry_answer()
print(result)