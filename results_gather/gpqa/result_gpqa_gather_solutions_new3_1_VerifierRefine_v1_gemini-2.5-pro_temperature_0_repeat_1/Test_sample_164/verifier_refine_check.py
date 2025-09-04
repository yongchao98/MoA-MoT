import re

def check_answer(question: str, final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the given chemistry question.

    The function evaluates the provided answer against established chemical principles
    for ethylene polymerization and tandem catalysis.

    Args:
        question: The question string (not used in this specific checker, but good practice).
        final_answer_text: The text containing the final answer, e.g., "<<<A>>>".

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """

    # Parse the letter from the final answer string, e.g., "<<<A>>>" -> "A"
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D. Got: {final_answer_text}"
    
    proposed_answer = match.group(1)

    # Ground truth analysis of the statements:
    # The core problem is creating "regular branches" from only ethylene. This requires
    # selective oligomerization (e.g., trimerization to 1-hexene) followed by copolymerization.

    # Statement A: "One can use a catalyst of a group VIa transition metal..."
    # Correctness: TRUE. Chromium (Cr) is in Group VIa (Group 6) and is the basis for the
    # most important industrial catalysts for selective ethylene trimerization to 1-hexene.
    # This is the key technology for the desired transformation.

    # Statement B: "Certain noble metal catalysts can be used but are too expensive."
    # Correctness: MISLEADING/LESS ACCURATE. While noble metals (e.g., Pd) can make branched
    # polyethylene via "chain-walking", this does not typically produce "regular" branches of a
    # uniform length. It's not the primary method for the specific goal.

    # Statement C: "Such combined systems are already implemented on an industrial scale in the US."
    # Correctness: AMBIGUOUS/INCORRECT. The term "combined system" is key. While the component
    # technologies (oligomerization and polymerization) are industrialized in the US, they are
    # typically run in separate processes/plants. A true single-reactor tandem system is not
    # the standard widespread industrial practice. So, the statement is incorrect under a strict
    # technical interpretation.

    # Statement D: "Aluminum-based activators do not work for the essential additional reaction step."
    # Correctness: FALSE. The "essential additional reaction step" is oligomerization. The catalysts
    # for this (e.g., Cr-based) are almost universally activated by aluminum-based compounds
    # like MAO, MMAO, or trialkylaluminums. They are essential, not non-functional.

    # The single best answer is A, as it describes the core, unambiguous chemical principle.
    correct_answer = 'A'

    if proposed_answer == correct_answer:
        return "Correct"
    else:
        reasons = {
            'B': "The proposed answer B is incorrect. While the statement is plausible, noble metal catalysts typically produce hyperbranched structures via 'chain-walking', not the 'regular branches' specified in the question. The most accurate answer, A, describes the primary industrial technology for creating regular branches.",
            'C': "The proposed answer C is incorrect. The statement is ambiguous and likely false under a strict interpretation. Industrial processes typically separate the oligomerization and polymerization steps into different units, so a true 'combined system' in a single reactor is not the standard industrial implementation.",
            'D': "The proposed answer D is incorrect. The statement is factually false. The essential additional reaction step (ethylene oligomerization) requires aluminum-based activators (like MAO) to function; they are crucial, not non-functional."
        }
        # Provide a generic error if the proposed answer is A but the logic determined another answer was correct
        # (This case shouldn't be hit given the analysis, but it's good practice)
        if proposed_answer not in reasons:
             return f"The proposed answer {proposed_answer} is incorrect. The correct answer is {correct_answer} because it describes the fundamental and correct chemical principle for the desired synthesis."

        return reasons[proposed_answer]

# Example usage with the provided data:
question_text = """
In an industrial research lab, a scientist performs ethylene polymerization with a homogeneous organometallic catalyst system, generating a polymer of high density. He intends to add a second catalyst system to introduce regular branches in the polymer backbone, also only using ethylene as the reactant. He consults a senior scientist, who gives the following statements. “Such combined systems are already implemented on an industrial scale in the US. One can use a catalyst of a group VIa transition metal in combination with specific activators. Aluminum-based activators do not work for the essential additional reaction step. Certain noble metal catalysts can be used but are too expensive.”
Which of these four statements is correct regarding the formation of a polymer with regular branches using only ethylene as the monomer and a dual catalyst system?

A) One can use a catalyst of a group VIa transition metal in combination with specific activators.
B) Certain noble metal catalysts can be used but are too expensive.
C) Such combined systems are already implemented on an industrial scale in the US.
D) Aluminum-based activators do not work for the essential additional reaction step.
"""

final_answer = """
Here is a step-by-step analysis of the question and the provided statements:

**Step 1: Deconstruct the Core Chemical Problem**
The scientist's goal is to produce a polymer with "regular branches" using only ethylene as the reactant. This requires a tandem catalysis process involving two distinct reactions:
1.  **Oligomerization:** A first catalyst (the "second catalyst system" to be added) must selectively convert some of the ethylene into a specific short-chain alpha-olefin. For example, it would trimerize three ethylene molecules into one molecule of 1-hexene. This selectivity is key to forming "regular" branches (e.g., all butyl branches from 1-hexene).
2.  **Copolymerization:** A second catalyst (the original one) must then copolymerize the remaining ethylene with the *in-situ* generated alpha-olefin to form the final branched polymer, which would be a type of Linear Low-Density Polyethylene (LLDPE).

**Step 2: Evaluate Each Statement**

*   **A) One can use a catalyst of a group VIa transition metal in combination with specific activators.**
    *   **Analysis:** Group VIa (or Group 6) of the periodic table contains Chromium (Cr). Chromium-based catalysts are the state-of-the-art industrial technology for the highly selective oligomerization of ethylene to 1-hexene (trimerization) and 1-octene (tetramerization). This is precisely the function required to generate the comonomer for regular branches. These catalysts indeed require specific activators, typically organoaluminum compounds.
    *   **Conclusion:** This statement is a direct, unambiguous, and factually correct description of the core chemistry required.

*   **B) Certain noble metal catalysts can be used but are too expensive.**
    *   **Analysis:** Noble metal catalysts (e.g., based on Palladium) can produce branched polyethylene from ethylene via a "chain-walking" mechanism. However, this process typically creates a hyperbranched structure with a variety of branch lengths (methyl, ethyl, etc.), which does not fit the description of "regular branches." Furthermore, their high cost is a major barrier to industrial use for commodity polymers.
    *   **Conclusion:** This statement is not the best description for the specific goal of creating *regular* branches and is thus a less relevant, distracting option.

*   **C) Such combined systems are already implemented on an industrial scale in the US.**
    *   **Analysis:** This statement is ambiguous. While the core technologies (selective ethylene oligomerization and ethylene polymerization) are certainly implemented on an industrial scale in the US by companies like Chevron Phillips and Sasol, they are typically run in separate, dedicated plants. A true "combined system" in a single reactor is a major research goal but is not the standard, widespread industrial practice.
    *   **Conclusion:** Due to the ambiguity of the term "combined system," this statement is potentially misleading or incorrect under a strict technical interpretation.

*   **D) Aluminum-based activators do not work for the essential additional reaction step.**
    *   **Analysis:** The "essential additional reaction step" is the oligomerization of ethylene. The chromium-based catalysts described in statement A are almost universally activated by aluminum-based compounds, such as methylaluminoxane (MAO) or trialkylaluminums. Far from "not working," these activators are essential for the catalyst's function.
    *   **Conclusion:** This statement is factually incorrect.

**Step 3: Determine the Best Answer**
- Statement D is definitively false.
- Statement B is not the best fit for the problem's specific goal ("regular branches").
- Statement C is factually ambiguous and likely incorrect under a strict interpretation of "combined system."
- Statement A is a direct, unambiguous, and fundamentally correct statement describing the key chemical technology that enables the desired synthesis. It is the most accurate and robust statement among the choices.

<<<A>>>
"""

# print(check_answer(question_text, final_answer))