import re

def check_exoplanet_density_answer(llm_response_text: str) -> str:
    """
    Checks the correctness of the final answer for the exoplanet density question.

    The function verifies the answer based on the physical principle that for planets
    of the same composition, higher mass leads to greater gravitational compression
    and thus higher density.

    Args:
        llm_response_text: The full text of the LLM's response, which includes
                           the final answer in the format <<<X>>>.

    Returns:
        "Correct" if the answer is correct.
        A string explaining the error if the answer is incorrect.
    """
    # Step 1: Define the problem's constraints and physical principles.
    # We can represent the densities qualitatively based on the descriptions.
    # Let's assign a rank to each planet's density. A higher rank means higher density.
    # - a) Earth-mass, Earth-radius: Baseline density.
    # - b) 2 Earth-masses, ~5.5 g/cm^3: Density is explicitly stated to be same as baseline.
    # - d) 0.5 Earth-mass, same composition: Less mass -> less compression -> lower density.
    # - c) 5 Earth-masses, same composition: More mass -> more compression -> higher density.
    
    density_ranking = {
        'a': 2,  # Baseline
        'b': 2,  # Stated to be same as baseline
        'c': 3,  # Highest density
        'd': 1   # Lowest density
    }

    # Step 2: Determine the correct planet label based on our physical model.
    # The planet with the highest density is the one with the highest rank.
    correct_planet_label = max(density_ranking, key=density_ranking.get)

    # Step 3: Define the mapping from the final answer choices (A, B, C, D) to the planet descriptions (a, b, c, d).
    answer_mapping = {
        'A': 'c',
        'B': 'a',
        'C': 'b',
        'D': 'd'
    }

    # Step 4: Extract the provided answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_response_text)
    if not match:
        return "Error: The final answer in the format <<<X>>> was not found in the provided text."

    provided_answer_letter = match.group(1)
    
    # Step 5: Check if the provided answer letter corresponds to the correct planet label.
    if answer_mapping.get(provided_answer_letter) == correct_planet_label:
        return "Correct"
    else:
        # Find the correct answer letter for the error message.
        correct_answer_letter = [k for k, v in answer_mapping.items() if v == correct_planet_label][0]
        
        reason = (
            f"The provided answer is '{provided_answer_letter}', which corresponds to planet description '{answer_mapping.get(provided_answer_letter)}'.\n"
            f"This is incorrect. The correct answer is planet '{correct_planet_label}', which has the highest density.\n"
            f"Reasoning: For planets of the same composition, density increases with mass due to gravitational compression. "
            f"Planet 'c' (5 Earth masses) is the most massive and therefore will be the densest.\n"
            f"The correct final answer choice is '{correct_answer_letter}'."
        )
        return reason

# The user-provided response to be checked.
llm_answer_to_check = """
Here is a step-by-step analysis of the question and the provided answers.

**Step 1: Understand the Core Scientific Principle**
The question is about planetary density. Density is mass divided by volume (ρ = M/V). For planets, especially rocky ones with similar compositions, density is not constant. A key factor is **gravitational compression**. As a planet's mass increases, its own gravity becomes stronger. This stronger gravity compresses the materials the planet is made of (like its iron core and silicate mantle), forcing them into a smaller volume than they would otherwise occupy. This effect leads to an increase in the planet's overall average density.

**Step 2: Analyze Each Option Based on this Principle**

*   **a) An Earth-mass and Earth-radius planet:** This is our baseline reference, Earth itself. Its average density is approximately 5.51 g/cm³.

*   **b) A planet with 2 Earth masses and a density of approximately 5.5 g/cm³:** The density is explicitly given as ~5.5 g/cm³. This is essentially the same as Earth's density.

*   **c) A planet with the same composition as Earth but 5 times more massive than Earth:** This planet is a "Super-Earth". Because it has the same composition as Earth but is much more massive, it will experience significantly more gravitational compression. This will squeeze its materials much more tightly, resulting in an overall density that is **higher** than Earth's.

*   **d) A planet with the same composition as Earth but half the mass of Earth:** This is a "sub-Earth". With less mass, its self-gravity is weaker. This results in less gravitational compression compared to Earth. Therefore, its overall density will be **lower** than Earth's.

**Step 3: Compare the Densities**
Based on the analysis:
*   Density of planet (d) < Density of Earth (~5.51 g/cm³)
*   Density of planet (a) ≈ 5.51 g/cm³
*   Density of planet (b) ≈ 5.5 g/cm³
*   Density of planet (c) > 5.51 g/cm³

The logical conclusion is that the planet described in option (c) has the highest density.

**Step 4: Select the Correct Final Answer**
The question asks to choose from the final options A, B, C, and D, which are mapped as follows:
A) c
B) a
C) b
D) d

Since our analysis determined that the planet in description **(c)** is the densest, the correct final answer is **A**.

**Careful Points and Analysis of Candidate Answers:**
*   Almost all candidate answers correctly identify the physical principle that higher mass leads to higher density due to gravitational compression (for the same composition).
*   The primary source of error among the candidate answers is a failure to correctly map the logical choice (c) to the final answer letter (A). For example, Answer 2 correctly reasons that (c) is the densest but incorrectly outputs <<<C>>>, which corresponds to option (b).
*   The correct analysis involves both understanding the physics and carefully reading the question's format to select the right letter.

<<<A>>>
"""

# Run the check
result = check_exoplanet_density_answer(llm_answer_to_check)
print(result)