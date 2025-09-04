import re

def check_correctness(question_text, final_answer_text):
    """
    Checks the correctness of the final answer for the given physics question.

    The question asks to identify the particle NOT associated with spontaneously-broken symmetry (SSB).
    This means identifying the particle that is NOT a Goldstone boson.

    - Magnon: Goldstone boson of broken rotational symmetry in ferromagnets.
    - Pion: Pseudo-Goldstone boson of broken chiral symmetry in QCD.
    - Phonon: Goldstone boson of broken translational symmetry in crystals.
    - Skyrmion: A topological soliton, NOT a Goldstone boson.

    Therefore, the correct answer is the Skyrmion.
    """

    # Define the physical properties of each particle in the context of the question.
    # The question asks for the particle NOT associated with SSB (i.e., not a Goldstone boson).
    particle_info = {
        "Magnon": {
            "is_goldstone_boson": True,
            "description": "a Goldstone boson of broken rotational symmetry."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "description": "a pseudo-Goldstone boson of broken chiral symmetry."
        },
        "Skyrmion": {
            "is_goldstone_boson": False,
            "description": "a topological soliton, not a Goldstone boson."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "description": "a Goldstone boson of broken translational symmetry."
        }
    }

    # Parse the options from the question text to handle different orderings.
    # Example line: "A) Magnon"
    options_regex = r"([A-D])\)\s*(\w+)"
    matches = re.findall(options_regex, question_text)
    
    if not matches or len(matches) != 4:
        return "Error: Could not parse the options A, B, C, D from the question text."
        
    question_options = {opt[0]: opt[1] for opt in matches}

    # Find the correct option letter based on our physics knowledge.
    correct_particle_name = "Skyrmion"
    correct_option_letter = None
    for letter, particle_name in question_options.items():
        if particle_name == correct_particle_name:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Error: The correct answer '{correct_particle_name}' was not found in the question options."

    # Extract the chosen answer from the final response.
    answer_match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not answer_match:
        return "Incorrect: The final answer is not in the required format '<<<X>>>'."

    chosen_option_letter = answer_match.group(1)

    # Compare the chosen answer with the correct answer.
    if chosen_option_letter == correct_option_letter:
        return "Correct"
    else:
        chosen_particle_name = question_options.get(chosen_option_letter, "Unknown")
        chosen_particle_info = particle_info.get(chosen_particle_name, {})
        
        if chosen_particle_info.get("is_goldstone_boson"):
            reason = (f"The chosen particle, '{chosen_particle_name}', is incorrect because it *is* associated with "
                      f"spontaneously-broken symmetry. It is {chosen_particle_info['description']} "
                      f"The question asks for the particle that is *not* associated with SSB.")
        else:
            # This case should not be reached if the logic is sound, but as a fallback.
            reason = (f"The chosen option '{chosen_option_letter}' is incorrect. "
                      f"The correct option is '{correct_option_letter}' ({correct_particle_name}).")

        return f"Incorrect: {reason}"

# The original question text to parse the options correctly.
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Which of the following (effective) particles is not associated with a spontaneously-broken symmetry?

A) Magnon
B) Pion
C) Skyrmion
D) Phonon

Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response, e.g. <<<C>>>, <<<A>>>
"""

# The final answer provided by the LLM.
llm_answer = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate answers and a review of the core physics concepts.

### Step 1: Understand the Core Physical Principle

The question asks to identify which particle is *not* associated with a spontaneously-broken symmetry (SSB). The most direct and fundamental association between a particle and SSB is described by **Goldstone's Theorem**. This theorem states that when a continuous global symmetry is spontaneously broken, the system must exhibit massless (or, if the symmetry is only approximate, low-mass) particles. These particles are called **Nambu-Goldstone bosons**.

Therefore, the question can be rephrased as: "Which of the following particles is not a Nambu-Goldstone boson?" All 15 candidate answers correctly adopt this interpretation.

### Step 2: Analyze Each Particle Based on Consensus and Fact-Checking

There is a complete consensus in the physical descriptions across all 15 candidate answers, which is confirmed by established physics principles.

*   **Magnon**: In a ferromagnet, the underlying physical laws are rotationally symmetric (no preferred direction for spin). However, the ground state of the material has all atomic spins aligned in a single direction, thus "spontaneously breaking" the rotational symmetry. The resulting low-energy excitations are spin waves, and their quanta are **magnons**. Magnons are the Goldstone bosons of this broken rotational symmetry. Thus, the magnon is directly associated with SSB.

*   **Pion**: In the theory of the strong force (Quantum Chromodynamics), there is an approximate symmetry called chiral symmetry. This symmetry is spontaneously broken by the vacuum. The resulting particles are the **pions**. Because the original symmetry is only approximate (it is explicitly broken by the small mass of quarks), the pions have a small mass and are technically called "pseudo-Goldstone bosons". Their existence is a classic example of a direct consequence of SSB.

*   **Phonon**: In a crystal, atoms arrange themselves into a fixed lattice. This breaks the continuous translational symmetry of free space (the laws of physics are the same everywhere, but the crystal itself is not). The collective vibrations of this lattice are quantized as **phonons**. The long-wavelength acoustic phonons are the Goldstone bosons corresponding to this broken translational symmetry. Thus, the phonon is directly associated with SSB.

*   **Skyrmion**: A skyrmion is a fundamentally different type of entity. It is a **topological soliton**, which is a stable, particle-like "knot" or "twist" in a field configuration. Its stability is guaranteed by its topology (a mathematical property related to its shape that cannot be smoothly "untied"), not by Goldstone's theorem. While skyrmions can arise in systems that also feature SSB (for instance, they can be modeled as configurations of pion fields), the skyrmion itself is not the Goldstone boson. It is a massive, complex, non-linear excitation, distinct from the massless (or light) Goldstone bosons.

### Step 3: Synthesize and Conclude

The analysis reveals a clear distinction:
*   **Associated with SSB (as Goldstone bosons):** Magnon, Pion, Phonon.
*   **Not associated with SSB (in the sense of being a Goldstone boson):** Skyrmion.

All three other particles are textbook examples of the direct particle manifestation of spontaneously broken continuous symmetries. The skyrmion is the outlier, belonging to a different class of excitations known as topological solitons.

### Step 4: Final Answer Selection

The question asks which particle is *not* associated with a spontaneously-broken symmetry. Based on the unanimous physical reasoning, this particle is the **Skyrmion**. According to the options provided in the question:

A) Magnon
B) Pion
C) Skyrmion
D) Phonon

The correct option is C.

<<<C>>>
"""

# Run the check
result = check_correctness(question, llm_answer)
print(result)