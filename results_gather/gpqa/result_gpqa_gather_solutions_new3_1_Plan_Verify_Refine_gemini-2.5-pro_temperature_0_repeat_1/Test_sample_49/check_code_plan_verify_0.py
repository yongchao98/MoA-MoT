import re

def check_correctness_of_feynman_loop_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the Feynman loop question.

    The core logic is based on the principle that in 4D Quantum Field Theory,
    each loop contributes a factor of 1/(4pi)^2 to the amplitude.
    The given expression has a factor of 1/(4pi)^6.
    Therefore, if L is the number of loops, (1/(4pi)^2)^L = 1/(4pi)^6.
    This implies 2L = 6, so L = 3.

    Args:
        llm_answer_text: The string containing the full response from the LLM,
                         including the final answer in <<<X>>> format.

    Returns:
        A string indicating "Correct" or the reason for being incorrect.
    """
    # --- 1. Define the problem's parameters and correct answer ---
    
    # Options given in the question
    options = {'A': 6, 'B': 2, 'C': 1, 'D': 3}
    
    # The key term in the expression is "1/(4pi)^6".
    # The total exponent of (4pi) in the denominator is 6.
    total_exponent = 6
    
    # In 4D QFT, each loop contributes a factor of 1/(4pi)^2.
    # The exponent contributed by each loop is 2.
    exponent_per_loop = 2
    
    # Calculate the correct number of loops
    try:
        correct_num_loops = total_exponent / exponent_per_loop
    except ZeroDivisionError:
        return "Error: Exponent per loop cannot be zero."

    # Find the correct option letter
    correct_option = None
    for opt, val in options.items():
        if val == correct_num_loops:
            correct_option = opt
            break
            
    if correct_option is None:
        return f"Internal check error: The calculated number of loops ({correct_num_loops}) does not match any of the options."

    # --- 2. Parse the LLM's provided answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    
    if not match:
        return "The answer is incorrect because it is not in the expected format '<<<X>>>' where X is one of A, B, C, or D."
    
    provided_option = match.group(1)
    provided_answer_value = options.get(provided_option)

    # --- 3. Compare and return the verdict ---
    if provided_option == correct_option:
        return "Correct"
    else:
        reason = (f"The answer is incorrect. "
                  f"The core principle is that each loop in 4D QFT contributes a factor of 1/(4pi)^2. "
                  f"The expression has a factor of 1/(4pi)^6. "
                  f"Therefore, the number of loops L is found by solving 2*L = 6, which gives L = 3. "
                  f"The correct answer is 3, which corresponds to option {correct_option}. "
                  f"The provided answer was option {provided_option}, which corresponds to {provided_answer_value} loops.")
        return reason

# The user-provided answer from the LLM
llm_answer = """
Here is a thinking process to arrive at the solution:

**Initial Plan:**
1.  Analyze the provided candidate answers to identify the core reasoning used.
2.  Verify the physical principle behind the reasoning.
3.  Apply the principle to the problem's expression to calculate the number of loops.
4.  Compare the result with the given options (A, B, C, D) and the candidates' final answers.
5.  Synthesize the analysis to produce the final, correct answer.

**Execution of Step 1 & 2 (Analysis and Verification):**
I will start by analyzing the reasoning provided in the candidate answers. A clear pattern emerges: almost all agents that provide a detailed explanation identify the `1/(4pi)^6` term as the key. They state that in four-dimensional spacetime, each loop (`L`) in a Feynman diagram contributes a factor of `1/(4pi)^2` to the final expression. This is a standard result from performing loop integrals using dimensional regularization. I can proceed with high confidence in this principle.

**Execution of Step 3 (Calculation):**
The next logical step is to perform the calculation based on this principle. I will use a simple Python script to formalize this.

*   **Principle:** Total factor = (Factor per loop) ^ (Number of loops)
*   **In this problem:** `1/(4pi)^6 = (1/(4pi)^2)^L`
*   **Equation for exponents:** `6 = 2 * L`
*   **Solve for L:** `L = 6 / 2 = 3`

I will now execute a Python script to confirm this and check it against the options.
