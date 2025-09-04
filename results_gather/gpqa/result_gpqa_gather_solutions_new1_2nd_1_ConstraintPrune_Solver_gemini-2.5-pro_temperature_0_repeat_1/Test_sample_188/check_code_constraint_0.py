import re

def check_answer_correctness(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given physics question.

    The question asks to identify the particle not associated with a spontaneously-broken
    symmetry, which in the context of Goldstone's theorem, means identifying the particle
    that is not a Goldstone boson.

    Args:
        question_text (str): The original question text with options.
        llm_answer_text (str): The full text of the LLM's response, including the final answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Step 1: Define the ground truth based on physics principles.
    # Goldstone's theorem states that spontaneously breaking a continuous symmetry
    # results in a massless particle called a Goldstone boson.
    # The question asks which particle is NOT such a boson.
    particle_info = {
        "Skyrmion": {
            "is_goldstone_boson": False,
            "type": "Topological Soliton",
            "description": "A stable, particle-like 'knot' in a field, whose stability is guaranteed by topology, not Goldstone's theorem."
        },
        "Phonon": {
            "is_goldstone_boson": True,
            "type": "Goldstone Boson",
            "description": "Results from the spontaneous breaking of continuous translational symmetry in a crystal."
        },
        "Magnon": {
            "is_goldstone_boson": True,
            "type": "Goldstone Boson",
            "description": "Results from the spontaneous breaking of continuous rotational symmetry in a ferromagnet."
        },
        "Pion": {
            "is_goldstone_boson": True,
            "type": "Pseudo-Goldstone Boson",
            "description": "Results from the spontaneous breaking of approximate chiral symmetry in QCD."
        }
    }

    # Step 2: Parse the options from the question text provided in the final answer block.
    # The final answer block provides the option mapping it used.
    options_regex = r"A\) (Skyrmion)\nB\) (Phonon)\nC\) (Magnon)\nD\) (Pion)"
    match = re.search(options_regex, llm_answer_text, re.IGNORECASE)
    
    if not match:
        return "Error: Could not parse the options list from the provided answer text. The format seems to have changed."

    options_map = {
        "A": match.group(1).strip(),
        "B": match.group(2).strip(),
        "C": match.group(3).strip(),
        "D": match.group(4).strip()
    }

    # Step 3: Determine the correct answer based on the ground truth.
    correct_particle_name = None
    for particle, info in particle_info.items():
        if not info["is_goldstone_boson"]:
            correct_particle_name = particle
            break
    
    if correct_particle_name is None:
        return "Error in checker logic: Could not determine the correct particle."

    correct_option_letter = None
    for letter, name in options_map.items():
        if name.lower() == correct_particle_name.lower():
            correct_option_letter = letter
            break

    if correct_option_letter is None:
        return f"Error in checker logic: The correct particle '{correct_particle_name}' was not found in the parsed options."

    # Step 4: Extract the LLM's final answer.
    llm_answer_match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not llm_answer_match:
        return "Error: Could not find the final answer in the required '<<<A>>>' format in the LLM's response."
    
    llm_chosen_option = llm_answer_match.group(1)

    # Step 5: Compare the LLM's answer with the correct answer and provide a verdict.
    if llm_chosen_option == correct_option_letter:
        return "Correct"
    else:
        llm_chosen_particle = options_map.get(llm_chosen_option, "Unknown")
        correct_particle_info = particle_info[correct_particle_name]
        llm_chosen_particle_info = particle_info.get(llm_chosen_particle, {})
        
        reason = (
            f"Incorrect. The LLM chose '{llm_chosen_option}' which corresponds to '{llm_chosen_particle}'.\n"
            f"The correct answer is '{correct_option_letter}' which corresponds to '{correct_particle_name}'.\n"
            f"Reasoning: The question asks for the particle NOT associated with a spontaneously-broken symmetry in the sense of being a Goldstone boson.\n"
            f"- {correct_particle_name} is the correct answer because it is a {correct_particle_info['type']}. {correct_particle_info['description']}\n"
            f"- Phonons, Magnons, and Pions are all Goldstone or pseudo-Goldstone bosons, which are direct consequences of SSB.\n"
            f"The LLM's choice, '{llm_chosen_particle}', is incorrect because it is a {llm_chosen_particle_info.get('type', 'N/A')} and is directly associated with SSB."
        )
        return reason

# The user's provided text for the final answer.
# Note: The question text is embedded within the final answer block.
final_answer_text = """
Based on a comprehensive analysis of the provided candidate answers, here is a step-by-step breakdown to determine the final, correct answer.

### Step 1: Understand the Core Physical Principle

The question asks to identify which particle is *not* associated with a spontaneously-broken symmetry (SSB). The key physical principle connecting particles to SSB is **Goldstone's Theorem**. This theorem states that for every continuous global symmetry that is spontaneously broken, a massless (or very low-mass) particle, called a **Goldstone boson**, must emerge.

Therefore, the question can be rephrased as: "Which of the following particles is not a Goldstone boson (or pseudo-Goldstone boson)?"

### Step 2: Analyze the Consensus Reasoning from Candidate Answers

All 15 candidate answers, despite having different option orderings and occasional errors in their final letter choice, are in **unanimous agreement** on the physical nature of each particle:

*   **Phonon:** Consistently identified as the Goldstone boson resulting from the spontaneous breaking of **continuous translational symmetry** when a crystal lattice forms from free space.
*   **Magnon:** Consistently identified as the Goldstone boson resulting from the spontaneous breaking of **continuous rotational symmetry** when spins align in a ferromagnet.
*   **Pion:** Consistently identified as a **pseudo-Goldstone boson** resulting from the spontaneous breaking of an approximate **chiral symmetry** in Quantum Chromodynamics (QCD).
*   **Skyrmion:** Consistently identified as a **topological soliton**, a fundamentally different class of object from a Goldstone boson.

### Step 3: Identify Careful Points and Nuances

1.  **The Meaning of "Associated With":** The most crucial point is the distinction between being *the result of* SSB versus simply *existing in a system with* SSB. Phonons, magnons, and pions are the direct particle consequences of SSB as described by Goldstone's theorem. A skyrmion, while it can be formed from fields of Goldstone bosons (like pions), is not the Goldstone boson itself. It is a more complex, massive, non-linear excitation.

2.  **Pion as a "Pseudo"-Goldstone Boson:** The pion is not strictly massless. This is because the chiral symmetry it originates from is only approximate—it is also explicitly broken by the small masses of quarks. This makes the pion a *pseudo*-Goldstone boson. However, its existence and relatively low mass are still a direct consequence of SSB, so it is correctly categorized with the others.

3.  **Skyrmion as a Topological Soliton:** This is the key differentiator. A skyrmion is a stable, particle-like "knot" or "twist" in a field. Its stability is guaranteed by a topological property (its "winding number"), not by Goldstone's theorem. It is a non-perturbative object, unlike Goldstone bosons which describe the low-energy, perturbative fluctuations of the system.

### Step 4: Synthesize and Conclude

The unanimous physical reasoning from all candidate answers leads to a clear conclusion: Phonons, Magnons, and Pions are all examples of Goldstone or pseudo-Goldstone bosons. Their existence is a direct manifestation of spontaneously broken continuous symmetries.

The **Skyrmion** is the outlier. It is a topological soliton, not a Goldstone boson. Therefore, it is the particle that is not associated with a spontaneously-broken symmetry in the sense of being its direct particle consequence.

To determine the final answer, we refer to the option list provided in the user's original question:
A) Skyrmion
B) Phonon
C) Magnon
D) Pion

The correct particle is the Skyrmion, which corresponds to option A.

<<<A>>>
"""

# The question text is not strictly needed as the final answer block restates the options,
# but we pass it for completeness.
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Which of the following (effective) particles is not associated with a spontaneously-broken symmetry?

A) Skyrmion
B) Phonon
C) Magnon
D) Pion
"""

# Run the check
result = check_answer_correctness(question, final_answer_text)
print(result)