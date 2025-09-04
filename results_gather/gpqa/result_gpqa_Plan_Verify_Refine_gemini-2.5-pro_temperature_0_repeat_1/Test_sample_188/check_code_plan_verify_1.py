import re

def check_particle_physics_answer(llm_answer_string: str) -> str:
    """
    Checks the answer to a physics question about spontaneous symmetry breaking.

    The function uses a knowledge base to determine the correct answer and compares
    it with the provided answer from an LLM.

    Question: Which of the following (effective) particles is not associated 
              with a spontaneously-broken symmetry?
              A) Pion
              B) Skyrmion
              C) Magnon
              D) Phonon
    """

    # Knowledge base defining the relationship of each particle to Spontaneous Symmetry Breaking (SSB).
    # The key distinction is whether the particle is a Goldstone boson, the direct quantum
    # associated with a broken continuous symmetry.
    particle_knowledge = {
        'A': {
            'name': 'Pion',
            'is_goldstone_boson': True,
            'explanation': 'Pions are the (pseudo-)Goldstone bosons of spontaneously broken chiral symmetry in QCD. Thus, they are directly associated with SSB.'
        },
        'B': {
            'name': 'Skyrmion',
            'is_goldstone_boson': False,
            'explanation': 'A Skyrmion is a topological soliton, not a Goldstone boson. While the theory describing Skyrmions (the Skyrme model) has a spontaneously broken symmetry, the Skyrmion itself is a massive, stable configuration of the fields, distinct from the massless Goldstone boson excitations (pions).'
        },
        'C': {
            'name': 'Magnon',
            'is_goldstone_boson': True,
            'explanation': 'Magnons are the Goldstone bosons that arise from the spontaneous breaking of global spin-rotation symmetry in ferromagnets. Thus, they are directly associated with SSB.'
        },
        'D': {
            'name': 'Phonon',
            'is_goldstone_boson': True,
            'explanation': 'Phonons are the Goldstone bosons that arise from the spontaneous breaking of continuous translational symmetry by a crystal lattice. Thus, they are directly associated with SSB.'
        }
    }

    # The question asks which particle is "not associated with" SSB. In this context,
    # the most precise interpretation is to identify the particle that is NOT a Goldstone boson.
    # Pions, Magnons, and Phonons are all Goldstone bosons. The Skyrmion is not.
    correct_option = 'B'

    # Extract the letter from the LLM's answer using regex.
    # The format is assumed to be <<<X>>> where X is the letter choice.
    match = re.search(r'<<<([A-D])>>>', llm_answer_string)
    
    if not match:
        return "Error: The answer is not in the required format. Please provide the answer as <<<A>>>, <<<B>>>, <<<C>>>, or <<<D>>>."

    user_choice = match.group(1)

    if user_choice == correct_option:
        return "Correct"
    else:
        # If the choice is incorrect, provide a detailed explanation.
        wrong_choice_info = particle_knowledge[user_choice]
        correct_choice_info = particle_knowledge[correct_option]
        
        reason = (
            f"Incorrect. The provided answer was <<<{user_choice}>>>.\n"
            f"Reason: The chosen particle, the {wrong_choice_info['name']}, is incorrect because it IS a direct consequence (a Goldstone boson) of spontaneous symmetry breaking. "
            f"{wrong_choice_info['explanation']}\n"
            f"The question asks for the particle that is NOT associated in this direct way.\n"
            f"The correct answer is B) Skyrmion. {correct_choice_info['explanation']}"
        )
        return reason

# Example usage:
# Let's assume the LLM provided the following response text, which includes the final answer.
# Note: The provided text in the prompt is just a justification and does not contain the final answer.
# We will test the checker with a hypothetical correct answer.
llm_response_correct = """
Based on the analysis, the pion, magnon, and phonon are all Goldstone bosons resulting from spontaneously broken symmetries. The skyrmion, however, is a topological soliton, a different class of object. Therefore, it is the correct answer.
<<<B>>>
"""
print(check_particle_physics_answer(llm_response_correct))

# Now, let's test with a hypothetical incorrect answer.
llm_response_incorrect = """
Pions are pseudo-Goldstone bosons, meaning the symmetry is only approximate. Perhaps this means they are not truly associated with it.
<<<A>>>
"""
print(check_particle_physics_answer(llm_response_incorrect))