import re

def check_physics_question_answer(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer to a physics question about spontaneous symmetry breaking.

    The function uses a knowledge base to determine the relationship of each particle
    to Spontaneous Symmetry Breaking (SSB) and Goldstone's Theorem.
    """

    # Knowledge base defining the nature of each particle in relation to SSB.
    # The primary distinction is whether a particle is a Goldstone boson, which is the
    # direct consequence of a spontaneously broken continuous symmetry.
    particle_knowledge = {
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "A topological soliton. It is a stable, non-perturbative configuration of a field, not the massless excitation (Goldstone boson) that arises from SSB."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "description": "A pseudo-Goldstone boson. It results from the spontaneous breaking of chiral symmetry in QCD. It would be massless if not for the small, explicit symmetry breaking by quark masses."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson resulting from the spontaneous breaking of continuous translational symmetry when a liquid crystallizes into a solid."
        },
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson resulting from the spontaneous breaking of spin-rotational symmetry in a ferromagnet."
        }
    }

    options_map = {"A": "Skyrmion", "B": "Pion", "C": "Phonon", "D": "Magnon"}

    # Determine the correct answer based on the knowledge base.
    # The question asks which particle is NOT associated with SSB. In this context,
    # it means which particle is not a Goldstone boson.
    correct_particle_name = None
    for particle, data in particle_knowledge.items():
        if not data["is_goldstone_boson"]:
            correct_particle_name = particle
            break
    
    if not correct_particle_name:
        return "Error in checker logic: Could not identify the correct answer from the knowledge base."

    correct_option = [opt for opt, name in options_map.items() if name == correct_particle_name][0]

    # Extract the LLM's answer from the response string.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Constraint not satisfied: The final answer was not provided in the required format, e.g., <<<A>>>."

    llm_option = match.group(1)
    
    # Compare the LLM's answer with the correct answer.
    if llm_option == correct_option:
        return "Correct"
    else:
        llm_particle_name = options_map.get(llm_option, "Unknown")
        reason = (
            f"The answer '{llm_option}' is incorrect. The particle '{llm_particle_name}' is indeed a direct result of spontaneous symmetry breaking.\n"
            f"Explanation: {particle_knowledge[llm_particle_name]['description']}\n"
            f"The correct answer is '{correct_option}) {correct_particle_name}', because it is a {particle_knowledge[correct_particle_name]['description'].lower()}"
        )
        return reason

# The provided response from the LLM to be checked
llm_response_to_check = """
**Proposer**: My initial hypothesis is that Pions, Phonons, and Magnons are all Goldstone bosons that arise from spontaneously broken continuous symmetries, while Skyrmions are topological solitons and thus are not. I will perform a search to confirm the nature of each particle and its direct association with spontaneous symmetry breaking.

**Skeptic**: This seems plausible, but we need to be precise. Are Pions true Goldstone bosons, or are they pseudo-Goldstone bosons, and does that distinction matter for this question? Is it possible for a topological object like a Skyrmion to also be considered a consequence of spontaneous symmetry breaking in some context? We must verify the defining characteristic of each particle in relation to this specific mechanism.

**Action**: I will execute a search to determine the relationship of each of the four particles (Pion, Phonon, Magnon, Skyrmion) to the phenomenon of spontaneously broken symmetry.
Spontaneous symmetry breaking occurs when the underlying laws of a system are symmetric, but the system's ground state is not. According to Goldstone's theorem, when a continuous symmetry is spontaneously broken, a massless particle, known as a Goldstone boson, must exist.

Here is an analysis of each particle:

*   **Pion (B):** Pions are the pseudo-Goldstone bosons that arise from the spontaneous breaking of the chiral-flavor symmetry of Quantum Chromodynamics (QCD). They would be massless if the quarks were massless, but because quarks have a small mass which explicitly breaks the symmetry, pions acquire a small mass and are therefore known as "pseudo-Goldstone bosons".
*   **Phonon (C):** Phonons are the Goldstone bosons that result from the spontaneous breaking of translational and rotational symmetry when a material crystallizes from a liquid to a solid state.
*   **Magnon (D):** Magnons, or spin waves, are the Goldstone bosons associated with the spontaneous breaking of rotational symmetry in a ferromagnet, where all spins align in a particular direction.
*   **Skyrmion (A):** A skyrmion is a topologically non-trivial, particle-like excitation in a field. While the existence of skyrmions is often due to underlying broken symmetries (for example, broken inversion symmetry can lead to interactions that stabilize skyrmions), the skyrmion itself is a topological soliton, not a Goldstone boson. A Goldstone boson is the massless (or light) excitation *resulting* from the symmetry breaking, whereas a skyrmion is a stable, localized topological configuration.

Therefore, the pion, phonon, and magnon are all directly associated with spontaneously-broken symmetries as the resulting Goldstone (or pseudo-Goldstone) bosons. The skyrmion is not.

<<<A>>>
"""

# Run the checker
result = check_physics_question_answer(llm_response_to_check)
print(result)