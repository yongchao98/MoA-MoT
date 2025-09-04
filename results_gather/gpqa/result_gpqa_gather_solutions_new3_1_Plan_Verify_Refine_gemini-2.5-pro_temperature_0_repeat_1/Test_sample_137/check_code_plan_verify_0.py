import re

def check_correctness(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given chemical reaction question.

    The function validates the answer based on fundamental principles of chemical kinetics:
    1.  **Arrhenius Equation/Collision Theory**: Reaction rates generally increase with temperature.
    2.  **Le Chatelier's Principle/Catalysis**: The concentration of reactants and catalysts directly affects the reaction rate. A decrease in catalyst concentration slows the reaction.
    3.  **pH Scale**: pH is a logarithmic measure of H+ concentration. An increase in pH from 1 to 4 means a 1000-fold decrease in [H+].
    4.  **Pressure Effects**: Pressure primarily affects reactions involving gases.

    Args:
        question (str): The question text.
        llm_answer (str): The LLM's answer text, including reasoning and the final choice.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """

    # --- Step 1: Define the ground truth based on chemical principles ---
    correct_option_letter = 'C'
    options = {
        'A': "The increased temperature of the solution",
        'B': "The increased volume of the solution",
        'C': "The increased pH of the solution",
        'D': "The increased pressure of the solution"
    }

    # --- Step 2: Parse the LLM's answer ---
    try:
        match = re.search(r'<<<([A-D])>>>', llm_answer)
        if not match:
            return "The answer is in an invalid format. The final choice in the format <<<X>>> is missing."
        
        selected_option_letter = match.group(1)
        reasoning_text = llm_answer.lower()
    except Exception as e:
        return f"Error parsing the LLM's answer: {e}"

    # --- Step 3: Evaluate the LLM's answer against the ground truth ---

    # Constraint 1: The selected option must be correct.
    if selected_option_letter != correct_option_letter:
        error_message = (
            f"The answer is incorrect. The selected option is {selected_option_letter} ('{options[selected_option_letter]}'), "
            f"but the correct option is {correct_option_letter} ('{options[correct_option_letter]}').\n"
            f"Reasoning: The problem states the reaction rate *slowed down*. An increase in temperature (Option A) would speed up the reaction. "
            f"A change in pressure (Option D) is irrelevant for this type of reaction. "
            f"The most significant factor is the increase in pH (Option C), which represents a 1000-fold decrease in the H+ ion catalyst concentration, causing the reaction to slow down dramatically."
        )
        return error_message

    # Constraint 2: The reasoning must correctly explain why increased temperature is the wrong cause.
    # It must state that higher temperature increases the rate, which contradicts the observation.
    if not ("temperature" in reasoning_text and ("increase" in reasoning_text or "speed up" in reasoning_text) and "rate" in reasoning_text and ("contradict" in reasoning_text or "opposition to" in reasoning_text or "cannot be the cause" in reasoning_text)):
        return "The answer is incorrect. Although the correct option 'C' was chosen, the reasoning fails to correctly explain why increased temperature (Option A) is the wrong answer. A correct explanation must state that higher temperature increases the reaction rate, which contradicts the observation that the rate slowed down."

    # Constraint 3: The reasoning must correctly explain why increased pH is the correct cause.
    # It must link increased pH to a decrease in H+ (catalyst) concentration, which slows the rate.
    if not ("ph" in reasoning_text and ("decrease" in reasoning_text or "slower" in reasoning_text) and "rate" in reasoning_text and ("h+" in reasoning_text or "catalyst" in reasoning_text or "hydrogen ion" in reasoning_text)):
        return "The answer is incorrect. Although the correct option 'C' was chosen, the reasoning fails to correctly explain why increased pH is the right answer. A correct explanation must link the increased pH to a decrease in H+ ion (catalyst) concentration, which in turn slows the reaction rate."

    # If all constraints are satisfied
    return "Correct"
