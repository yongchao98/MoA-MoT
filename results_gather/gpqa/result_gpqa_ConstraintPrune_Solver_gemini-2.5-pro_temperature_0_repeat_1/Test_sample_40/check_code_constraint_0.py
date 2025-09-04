import math

def check_correctness():
    """
    Checks the correctness of the given answer for the relativistic speed of light problem.
    """
    # The answer provided by the LLM is C.
    llm_answer_option = 'C'

    # Define test parameters.
    # n = refractive index of the medium (must be > 1 for a physical medium like glass).
    # v = velocity of the medium relative to the observer (must be 0 <= v < 1, as c=1).
    n = 1.5
    v = 0.5
    c = 1.0

    # --- Ground Truth Calculation ---
    # According to physics, the speed of light in the rest frame of the glass (S') is u_prime.
    # u_prime = c / n
    u_prime = c / n

    # The speed 'u' as seen by the observer in the laboratory frame (S) is found
    # using the Einstein velocity-addition formula: u = (u_prime + v) / (1 + u_prime * v / c^2)
    expected_u = (u_prime + v) / (1 + u_prime * v / (c**2))

    # --- Evaluate the LLM's Answer ---
    # Define the formulas for all given options.
    options = {
        "A": lambda n, v: (1 + n * v) / (n - v),
        "B": lambda n, v: (1 - n * v) / (n + v),
        "C": lambda n, v: (1 + n * v) / (n + v),
        "D": lambda n, v: 1.0
    }

    # Check if the selected option exists.
    if llm_answer_option not in options:
        return f"Invalid option '{llm_answer_option}' provided. Options are A, B, C, D."

    # Calculate the result using the formula from the chosen answer.
    answer_formula = options[llm_answer_option]
    calculated_u = answer_formula(n, v)

    # --- Verification ---
    # Compare the ground truth with the result from the answer's formula.
    # Use math.isclose() for safe floating-point comparison.
    if math.isclose(expected_u, calculated_u):
        return "Correct"
    else:
        # If incorrect, provide a clear reason.
        # The algebraic simplification of the velocity-addition formula is:
        # u = (c/n + v) / (1 + (c/n)*v/c^2) = (1/n + v) / (1 + v/n)
        #   = n*(1/n + v) / n*(1 + v/n) = (1 + n*v) / (n + v)
        # This simplified formula is option C.
        correct_formula_str = "(1 + n*v) / (n + v)"
        
        return (f"Incorrect. The answer is given as option {llm_answer_option}, which corresponds to the formula "
                f"that yields a value of {calculated_u:.6f} for n={n} and v={v}. "
                f"However, the correct speed, derived from the relativistic velocity-addition formula, is {expected_u:.6f}. "
                f"The correct algebraic formula is {correct_formula_str}, which corresponds to option C.")

# Execute the check and print the result.
result = check_correctness()
print(result)