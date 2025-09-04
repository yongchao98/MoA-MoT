import re

def check_correctness(question, final_answer):
    """
    Checks the correctness of the final answer for the given question about ethylene polymerization.

    The function encodes established chemical facts as a knowledge base and evaluates each statement.
    It then compares its conclusion with the provided final answer.

    Args:
        question (str): The full text of the question.
        final_answer (str): The final answer provided by the LLM, including the reasoning.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """

    # --- Knowledge Base ---
    # This section encodes the established facts of the chemical process.

    # Fact 1: The "essential additional reaction step" is the selective oligomerization of ethylene
    # into a specific alpha-olefin (e.g., 1-hexene) to create "regular branches".
    # Fact 2: Chromium (Cr), a Group VIa (Group 6) metal, is the key component in state-of-the-art
    # industrial catalysts for this selective oligomerization.
    # Fact 3: These chromium catalysts require an activator (co-catalyst) to function.
    # Fact 4: The most common and essential activators are aluminum-based compounds like
    # methylaluminoxane (MAO) or trialkylaluminums.
    # Fact 5: Noble metal catalysts (e.g., Palladium) can produce branched polyethylene from ethylene,
    # but typically via a "chain-walking" mechanism, which results in a variety of branch lengths
    # (hyperbranching), not the "regular branches" specified.
    # Fact 6: True single-reactor, dual-catalyst "combined systems" are a major research goal but are
    # not the widespread industrial standard. Industrial practice often separates the oligomerization
    # and polymerization steps into different reactors or plants.

    analysis = {
        'A': {
            'is_correct': False,
            'reason': 'The statement claims "Such combined systems are already implemented on an industrial scale in the US." This is misleading. While the core technologies (oligomerization and polymerization) are industrial, true single-reactor "combined systems" are not the widespread standard. Industrial processes typically separate these steps.'
        },
        'B': {
            'is_correct': False,
            'reason': 'The statement claims "Aluminum-based activators do not work for the essential additional reaction step." This is factually incorrect. Aluminum-based activators like MAO are essential for the function of the chromium oligomerization catalysts.'
        },
        'C': {
            'is_correct': True,
            'reason': 'The statement claims "One can use a catalyst of a group VIa transition metal in combination with specific activators." This is correct. Chromium (a Group VIa metal) catalysts are the key technology for selectively producing alpha-olefins from ethylene, which is the required step to create regular branches.'
        },
        'D': {
            'is_correct': False,
            'reason': 'The statement claims "Certain noble metal catalysts can be used but are too expensive." While plausible, this is not the best description for the specific goal. Noble metal catalysts typically produce hyperbranched structures with irregular branches, which does not match the "regular branches" specified in the question.'
        }
    }

    # Find the correct answer based on the analysis
    correct_options = [opt for opt, res in analysis.items() if res['is_correct']]

    # Extract the letter from the final answer (e.g., '<<<C>>>' -> 'C')
    match = re.search(r'<<<([A-D])>>>', final_answer)
    if not match:
        return f"Invalid answer format. Could not find '<<<X>>>' in the provided answer."

    llm_answer_option = match.group(1)

    # Check if the LLM's answer is in the list of correct options
    if llm_answer_option in correct_options:
        # The LLM's reasoning should also align with the facts.
        # Let's check if the LLM correctly identified statement B as false.
        if "Statement B is definitively false" in final_answer or "statement is incorrect" in final_answer.lower() and "aluminum-based activators" in final_answer.lower():
             return "Correct"
        else:
             return f"The final answer choice '{llm_answer_option}' is correct, but the reasoning provided is flawed or incomplete. It fails to correctly identify that statement B is definitively false."

    else:
        # The LLM's answer is incorrect.
        llm_reason = analysis[llm_answer_option]['reason']
        correct_reason = analysis[correct_options[0]]['reason']
        return (f"Incorrect. The final answer was '{llm_answer_option}', but this is wrong. {llm_reason}\n"
                f"The correct answer is '{correct_options[0]}'. {correct_reason}")


# The user's provided question and final answer
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
In an industrial research lab, a scientist performs ethylene polymerization with a homogeneous organometallic catalyst system, generating a polymer of high density. He intends to add a second catalyst system to introduce regular branches in the polymer backbone, also only using ethylene as the reactant. He consults a senior scientist, who gives the following statements. “Such combined systems are already implemented on an industrial scale in the US. One can use a catalyst of a group VIa transition metal in combination with specific activators. Aluminum-based activators do not work for the essential additional reaction step. Certain noble metal catalysts can be used but are too expensive.”
Which of these four statements is correct regarding the formation of a polymer with regular branches using only ethylene as the monomer and a dual catalyst system?

A) Such combined systems are already implemented on an industrial scale in the US.
B) Aluminum-based activators do not work for the essential additional reaction step.
C) One can use a catalyst of a group VIa transition metal in combination with specific activators.
D) Certain noble metal catalysts can be used but are too expensive.
"""

final_answer_text = """
Here is a step-by-step analysis of the question to determine the correct statement.

**Step 1: Deconstruct the Core Problem**
The scientist's goal is to produce polyethylene with "regular branches" using only ethylene as the reactant. This requires a tandem catalysis process within a dual-catalyst system.
*   **Catalyst 1 (Polymerization):** Takes ethylene and polymerizes it into long chains. This catalyst must also be able to incorporate alpha-olefins.
*   **Catalyst 2 (Oligomerization):** This is the "second catalyst system" to be added. Its function is to perform the "essential additional reaction step," which is the selective oligomerization of ethylene into a specific short-chain alpha-olefin (e.g., trimerizing three ethylene molecules into one 1-hexene molecule).
The 1-hexene then acts as a comonomer, which Catalyst 1 incorporates into the main polyethylene chain, creating regular butyl branches.

**Step 2: Evaluate Each Statement**

*   **A) Such combined systems are already implemented on an industrial scale in the US.**
    This statement's accuracy depends on the interpretation of "combined systems." While the core technology (selective oligomerization of ethylene to 1-hexene) is indeed industrialized in the US by companies like Chevron Phillips and Sasol, it is typically done in dedicated plants. The alpha-olefins are then used as comonomers in separate polymerization reactors. A true single-reactor, dual-catalyst system for this purpose is not the widespread industrial standard. Therefore, this statement is likely incorrect under a strict technical interpretation.

*   **B) Aluminum-based activators do not work for the essential additional reaction step.**
    The "essential additional reaction step" is the oligomerization of ethylene. The catalysts used for this step, particularly the highly selective chromium-based ones, require an activator (or cocatalyst) to function. The scientific and industrial literature overwhelmingly shows that aluminum-based compounds, such as methylaluminoxane (MAO), modified MAO (MMAO), and various trialkylaluminum compounds, are the most common and essential activators for these systems. Therefore, this statement is factually incorrect.

*   **C) One can use a catalyst of a group VIa transition metal in combination with specific activators.**
    Group VIa (or Group 6) of the periodic table contains Chromium (Cr). Chromium-based catalysts are the state-of-the-art technology for the highly selective oligomerization of ethylene to 1-hexene (trimerization) and 1-octene (tetramerization). This is precisely the function required to generate the comonomer for regular branching. These catalysts indeed require "specific activators," which, as noted above, are typically aluminum-based. This statement accurately describes the fundamental chemistry of the most successful industrial solution to the problem. This statement is correct.

*   **D) Certain noble metal catalysts can be used but are too expensive.**
    While catalysts based on noble metals like palladium can produce branched polyethylene from ethylene via a "chain-walking" mechanism, this process typically creates a variety of branch types and lengths (hyperbranching), not the "regular branches" specified in the problem. Furthermore, the high cost of noble metals makes them economically unviable for producing a commodity polymer like polyethylene. So, while the statement is plausible in a general sense, it is not the correct or relevant technology for the specific goal described.

**Step 3: Final Conclusion**
Statement B is definitively false. Statement A is incorrect under a strict interpretation of industrial practice. Statement D describes a different, less relevant technology. Statement C provides a direct, accurate, and fundamental description of the key chemical technology used to achieve the scientist's goal. It is the most accurate and relevant statement among the four choices.

<<<C>>>
"""

# Run the check
result = check_correctness(question_text, final_answer_text)
print(result)