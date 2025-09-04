import sympy

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the underlying formula
    against fundamental physical principles using symbolic mathematics.
    """
    # Define symbolic variables for refractive index (n) and velocity (v)
    n, v = sympy.symbols('n v', real=True, positive=True)
    
    # The speed of light in vacuum is given as c=1
    c = 1

    # The options from the question are represented as symbolic expressions
    options = {
        'A': (1 - n*v) / (n + v),
        'B': (1 + n*v) / (n + v),
        'C': sympy.Integer(1),
        'D': (1 + n*v) / (n - v)
    }

    # The final answer provided by the LLM analysis is 'B'
    llm_answer_choice = 'B'
    
    # Get the formula corresponding to the LLM's answer
    formula_to_check = options.get(llm_answer_choice)

    if formula_to_check is None:
        return f"Error: The answer choice '{llm_answer_choice}' is not a valid option."

    # --- Constraint 1: The Vacuum Case ---
    # If the medium is a vacuum, its refractive index n = 1.
    # The speed of light in a vacuum must be c=1 for any observer.
    # We substitute n=1 into the formula and check if it simplifies to 1.
    speed_in_vacuum = sympy.simplify(formula_to_check.subs(n, 1))
    
    if speed_in_vacuum != c:
        return (f"Incorrect. The formula for answer {llm_answer_choice}, which is {formula_to_check}, "
                f"fails the vacuum constraint. When n=1, the speed should be {c}, "
                f"but the formula gives {speed_in_vacuum}.")

    # --- Constraint 2: The Stationary Medium Case ---
    # If the medium is stationary, its velocity v = 0.
    # The speed of light in a stationary medium is c/n = 1/n.
    # We substitute v=0 into the formula and check if it simplifies to 1/n.
    speed_in_stationary_medium = sympy.simplify(formula_to_check.subs(v, 0))

    if speed_in_stationary_medium != c/n:
        return (f"Incorrect. The formula for answer {llm_answer_choice}, which is {formula_to_check}, "
                f"fails the stationary medium constraint. When v=0, the speed should be 1/n, "
                f"but the formula gives {speed_in_stationary_medium}.")
                
    # --- Final Check: Relativistic Speed Limit ---
    # The calculated speed u must not exceed the speed of light c=1 for any physical v (0 <= v < 1) and n > 1.
    # We check the sign of (u - c). It must be <= 0.
    # u - c = (1 + n*v)/(n+v) - 1 = (1 + n*v - n - v)/(n+v) = ((v-1)*(n-1))/(n+v)
    # For v < 1, (v-1) is negative. For n > 1, (n-1) is positive. For v,n > 0, (n+v) is positive.
    # So, (negative * positive) / positive = negative. The speed is always less than c.
    # This confirms the formula is relativistically sound.

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)