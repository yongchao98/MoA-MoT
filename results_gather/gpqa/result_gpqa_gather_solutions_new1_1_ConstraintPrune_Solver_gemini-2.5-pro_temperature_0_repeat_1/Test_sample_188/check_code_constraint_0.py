import re

def check_correctness_of_particle_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to a particle physics question.

    The question asks to identify which particle from a list is NOT associated with
    a spontaneously-broken symmetry (SSB). The core physical principle is that
    particles directly resulting from SSB are Goldstone bosons (or pseudo-Goldstone bosons).
    The particle that is NOT a Goldstone boson is the correct answer.

    Args:
        llm_answer_text: The full text of the LLM's response, including the final answer.

    Returns:
        A string indicating "Correct" or explaining the reason for the error.
    """

    # Step 1: Define the ground truth based on established physics principles.
    # This knowledge base maps each particle to its classification.
    # The key constraint is whether the particle is a Goldstone boson.
    knowledge_base = {
        "Pion": {
            "is_goldstone_boson": True,
            "description": "A pseudo-Goldstone boson from the spontaneous breaking of approximate chiral symmetry."
        },
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "A topological soliton, whose stability comes from topology, not Goldstone's theorem."
        },
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson from the spontaneous breaking of rotational symmetry in a ferromagnet."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "A Goldstone boson from the spontaneous breaking of continuous translational symmetry in a crystal."
        }
    }

    # Step 2: Define the mapping from the question's options to the particle names.
    # This is based on the provided LLM analysis text.
    options_map = {
        "A": "Pion",
        "B": "Skyrmion",
        "C": "Magnon",
        "D": "Phonon"
    }

    # Step 3: Parse the final answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not provided in the required format '<<<X>>>'."

    proposed_answer_key = match.group(1)
    proposed_particle_name = options_map.get(proposed_answer_key)

    if not proposed_particle_name:
        return f"Incorrect: The answer key '{proposed_answer_key}' is not a valid option."

    # Step 4: Verify the proposed answer against the ground truth.
    # The question asks which particle is NOT associated with SSB.
    # Therefore, the correct answer must be a particle that is NOT a Goldstone boson.

    # Find the correct particle based on our knowledge base.
    correct_particle = None
    for particle, properties in knowledge_base.items():
        if not properties["is_goldstone_boson"]:
            correct_particle = particle
            break
    
    # Find the key for the correct particle
    correct_key = next((key for key, value in options_map.items() if value == correct_particle), None)


    # Check if the proposed answer matches the correct one.
    if proposed_particle_name == correct_particle:
        # The answer is correct. Let's also verify that the reasoning is sound.
        # The reasoning should correctly identify all other particles as Goldstone bosons.
        llm_reasoning = llm_answer_text
        all_others_correctly_identified = True
        for key, particle in options_map.items():
            if key != proposed_answer_key: # Check the other particles
                if not knowledge_base[particle]["is_goldstone_boson"]:
                    # This would be a contradiction in the problem itself, but good to check.
                    all_others_correctly_identified = False
        
        if all_others_correctly_identified:
            return "Correct"
        else:
            return "Incorrect: The final answer letter is correct, but the problem premise as understood from the reasoning is flawed."

    else:
        # The proposed answer is incorrect.
        return (f"Incorrect: The provided answer is {proposed_answer_key} ({proposed_particle_name}). "
                f"However, the {proposed_particle_name} IS a Goldstone (or pseudo-Goldstone) boson and therefore IS associated with spontaneously-broken symmetry. "
                f"The correct answer is {correct_key} ({correct_particle}), which is a topological soliton, not a Goldstone boson.")

# The user's provided response text to be checked.
llm_response_to_check = """
Here is a step-by-step analysis to decide the final answer.

### Step 1: Define the Core Physical Principle and Constraint

The question asks which particle is *not* associated with a spontaneously-broken symmetry (SSB). The central physical principle connecting particles to SSB is **Goldstone's Theorem**.

*   **Constraint:** Goldstone's theorem states that when a continuous global symmetry is spontaneously broken, a massless (or low-mass) particle called a **Goldstone boson** must appear. These bosons are the "effective particles" directly associated with SSB.
*   **Task Model:** The task is to test each particle against the constraint: "Is this particle a Goldstone boson (or pseudo-Goldstone boson)?" The particle that fails this test is the correct answer.

### Step 2: Prune the Options by Applying the Constraint

We will evaluate each option to see if it can be pruned as being a Goldstone boson.

*   **A) Pion:** All candidate answers unanimously identify the pion as the **pseudo-Goldstone boson** resulting from the spontaneous breaking of approximate chiral symmetry in Quantum Chromodynamics (QCD).
    *   *Careful Point:* It is a "pseudo" Goldstone boson because the symmetry is not perfect (due to non-zero quark masses), which gives the pion a small mass. However, its existence is a direct consequence of SSB.
    *   **Result:** The Pion is a particle associated with SSB. This option is pruned.

*   **C) Magnon:** All candidate answers unanimously identify the magnon as the **Goldstone boson** resulting from the spontaneous breaking of rotational symmetry in a ferromagnet (where all spins align in one direction).
    *   **Result:** The Magnon is a particle associated with SSB. This option is pruned.

*   **D) Phonon:** All candidate answers unanimously identify the (acoustic) phonon as the **Goldstone boson** resulting from the spontaneous breaking of continuous translational symmetry when a crystal lattice forms.
    *   **Result:** The Phonon is a particle associated with SSB. This option is pruned.

### Step 3: Identify and Analyze the Remaining Candidate

After pruning options A, C, and D, only option B remains.

*   **B) Skyrmion:** All candidate answers unanimously identify the skyrmion as a **topological soliton**.
    *   *Careful Point:* A skyrmion is a stable, particle-like "lump" or "twist" in a field. Its stability comes from its topology (a "winding number"), not from Goldstone's theorem. While skyrmions can be made of Goldstone boson fields (e.g., pions) and exist in systems with SSB, the skyrmion itself is a massive, non-perturbative object and is fundamentally different from a Goldstone boson.
    *   **Result:** The Skyrmion is not a Goldstone boson. It is not the direct particle manifestation of SSB in the sense of Goldstone's theorem.

### Step 4: Final Conclusion

The analysis shows a clear distinction: Pions, Magnons, and Phonons are all examples of Goldstone (or pseudo-Goldstone) bosons, which are the defining particle consequence of spontaneously broken continuous symmetries. The Skyrmion is a different class of object, a topological soliton. Therefore, the Skyrmion is the particle that is not associated with a spontaneously-broken symmetry in the manner described by Goldstone's theorem.

Based on the question's mapping (A: Pion, B: Skyrmion, C: Magnon, D: Phonon), the correct choice is B.

<<<B>>>
"""

# Run the check and print the result.
result = check_correctness_of_particle_answer(llm_response_to_check)
print(result)