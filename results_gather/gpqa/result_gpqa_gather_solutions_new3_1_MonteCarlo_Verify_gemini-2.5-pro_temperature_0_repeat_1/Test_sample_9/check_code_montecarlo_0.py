import re

def check_correctness(question_text, llm_answers_text, final_answer_text):
    """
    Checks the correctness of the final answer based on the question's logic.

    The logic is as follows:
    1. Establish a baseline density for an Earth-like planet (a).
    2. Note the given density for planet (b).
    3. Apply the principle of gravitational compression:
       - Planet (c), with 5x Earth mass and same composition, will be denser than Earth.
       - Planet (d), with 0.5x Earth mass and same composition, will be less dense than Earth.
    4. Determine which planet is the densest based on this ranking.
    5. Check if the provided final answer corresponds to the correct planet.
    """

    # --- Data Extraction ---
    # The mapping from final answer choices (A, B, C, D) to planet descriptions (a, b, c, d)
    # This is parsed from the question text.
    # Example: "A) c\nB) b\nC) d\nD) a"
    try:
        mapping_lines = re.findall(r'([A-D])\)\s*([a-d])', question_text)
        if not mapping_lines:
            return "Error: Could not parse the answer choice mapping (e.g., 'A) c') from the question text."
        answer_mapping = {choice: desc for choice, desc in mapping_lines}
    except Exception as e:
        return f"Error parsing the question's answer mapping: {e}"

    # Extract the final answer letter from the provided text (e.g., "<<<A>>>")
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<X>>>' where X is A, B, C, or D."
    
    provided_answer_choice = match.group(1)

    # --- Logical Analysis ---
    # We don't need exact values, just the relative order.
    # Let's assign scores to represent relative density.
    densities = {}

    # a) Baseline (Earth)
    densities['a'] = 5.51  # Earth's actual average density

    # b) Given density
    densities['b'] = 5.50  # Approximately 5.5 g/cm^3

    # c) Same composition, 5x mass -> higher compression -> higher density
    # We know it's > 5.51, so we can assign a higher score.
    densities['c'] = densities['a'] + 2 # A value guaranteed to be the highest

    # d) Same composition, 0.5x mass -> lower compression -> lower density
    # We know it's < 5.51, so we can assign a lower score.
    densities['d'] = densities['a'] - 1 # A value guaranteed to be lower

    # Determine the description of the densest planet from our analysis
    densest_planet_description = max(densities, key=densities.get)

    # Find the correct final answer choice letter based on our analysis
    correct_choice = None
    for choice, desc in answer_mapping.items():
        if desc == densest_planet_description:
            correct_choice = choice
            break
    
    # --- Verification ---
    if provided_answer_choice == correct_choice:
        return "Correct"
    else:
        reason = (
            f"The provided answer is '{provided_answer_choice}', which corresponds to planet description '{answer_mapping.get(provided_answer_choice)}'.\n"
            f"However, the correct answer is '{correct_choice}', which corresponds to planet description '{densest_planet_description}'.\n"
            f"Reasoning: The key principle is gravitational compression. For planets with the same composition, a more massive planet is denser. "
            f"Planet (c) is 5 times more massive than Earth (a), so it will be significantly denser. "
            f"Planet (d) is less massive, so it will be less dense. "
            f"Planet (b) has a density explicitly stated to be similar to Earth's. "
            f"Therefore, the ranking of densities is (c) > (a) ≈ (b) > (d), making (c) the planet with the highest density."
        )
        return reason

# --- Inputs for the checker ---
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Among the following exoplanets, which one has the highest density?

a) An Earth-mass and Earth-radius planet.
b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm^3.
c) A planet with the same composition as Earth but 5 times more massive than Earth.
d) A planet with the same composition as Earth but half the mass of Earth.

A) c
B) b
C) d
D) a
"""

# We don't need the candidate answers for this logical check, only the final answer.
candidate_answers = "..." # Not needed for this checker

# The final answer provided by the LLM being checked.
final_answer = """
Here is a step-by-step analysis to determine the correct answer.

**Step 1: Understand the core physical principle.**
The question revolves around the concept of planetary density. Density is mass divided by volume (ρ = M/V). For planets made of the same materials (i.e., same composition), the density is not constant. A more massive planet exerts a stronger gravitational force on itself. This phenomenon, known as **gravitational compression**, squeezes the planet's material into a smaller volume than it would otherwise occupy, thus increasing its overall density.

**Step 2: Analyze each of the described planets.**

*   **a) An Earth-mass and Earth-radius planet:** This is our baseline reference, a planet identical to Earth. Its average density is approximately 5.51 g/cm³.
*   **b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm³:** The density is explicitly given as ~5.5 g/cm³. This is essentially the same as Earth's density.
*   **c) A planet with the same composition as Earth but 5 times more massive than Earth:** This planet is a "Super-Earth". Because it has the same composition as Earth but is much more massive, it will experience significantly stronger gravitational compression. This will make it denser than Earth. Therefore, its density will be > 5.51 g/cm³.
*   **d) A planet with the same composition as Earth but half the mass of Earth:** This planet is a "sub-Earth". With less mass, its self-gravity is weaker, leading to less gravitational compression compared to Earth. Therefore, its density will be < 5.51 g/cm³.

**Step 3: Compare the densities.**
Based on the analysis, we can rank the densities:
Density(d) < Density(a) ≈ Density(b) < Density(c)

The planet described in (c) has the highest density.

**Step 4: Select the correct final answer option.**
The question provides the following choices:
A) c
B) b
C) d
D) a

Our analysis concluded that the planet described in (c) is the one with the highest density. This corresponds to option A.

**Careful Points and Synthesis of Candidate Answers:**
*   **Consensus on Physics:** Nearly all candidate answers correctly identify that for a planet of a given composition, density increases with mass due to gravitational compression.
*   **Consensus on Planet:** All reasoned answers correctly conclude that description (c) represents the densest planet.
*   **Source of Error in Candidates:** The primary discrepancy among the candidate answers is not in the scientific reasoning but in mapping the correct description (c) to the final letter choice (A, B, C, or D). Many candidates seem to have assumed a standard A=a, B=b, C=c mapping, which is incorrect for this specific question's format. The correct mapping is explicitly provided.

By following the correct physical principles and carefully matching the conclusion to the provided options, the answer is unequivocally A.

<<<A>>>
"""

# Run the check
result = check_correctness(question, candidate_answers, final_answer)
print(result)