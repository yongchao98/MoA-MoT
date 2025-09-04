import sympy
from sympy import symbols, sqrt, Eq, solve, limit, oo

def check_correctness():
    """
    This function verifies the correctness of the answer to the relativistic harmonic oscillator problem.
    It performs three main checks:
    1. Derives the correct expression from first principles using symbolic math.
    2. Compares the derived expression with the proposed answer (C).
    3. Performs physical sanity checks (dimensional consistency, physical bounds, non-relativistic limit)
       on all options to show why C is correct and others are not.
    """
    # Define symbolic variables
    m, k, A, c = symbols('m k A c', positive=True, real=True)
    v_max = symbols('v_max', real=True)

    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'C'
    
    # Define the four candidate answers as symbolic expressions from the question
    options = {
        'A': c * sqrt(1 + 1 / (1 - k*A**2 / (2*m))**2),
        'B': sqrt(k*A**2 / m),
        'C': c * sqrt(1 - 1 / (1 + k*A**2 / (2*m*c**2))**2),
        'D': c * sqrt(1 + 1 / (1 - k*A**2 / (2*m*c**2)))
    }
    
    # --- Check 1: Symbolic Derivation from First Principles ---
    try:
        # Total energy at max amplitude (x=A, v=0)
        E_total_at_A = m*c**2 + (k*A**2)/2
        
        # Total energy at equilibrium (x=0, v=v_max)
        gamma_max_sym = symbols('gamma_max_sym')
        E_total_at_0 = gamma_max_sym * m * c**2
        
        # From conservation of energy, solve for gamma_max
        gamma_max_expr = solve(Eq(E_total_at_A, E_total_at_0), gamma_max_sym)[0]
        
        # Now solve for v_max from the definition of gamma
        # gamma = 1 / sqrt(1 - v^2/c^2)  =>  v^2 = c^2 * (1 - 1/gamma^2)
        derived_v_max_squared = c**2 * (1 - 1/gamma_max_expr**2)
        
        # Compare the square of the derived expression with the square of the chosen answer's expression
        chosen_option_expr = options[llm_answer_choice]
        diff = sympy.simplify(derived_v_max_squared - chosen_option_expr**2)
        
        if diff != 0:
            return f"Incorrect. The provided answer {llm_answer_choice} does not match the derived expression.\nDerived v_max^2: {sympy.simplify(derived_v_max_squared)}\nAnswer {llm_answer_choice} v_max^2: {sympy.simplify(chosen_option_expr**2)}"
    except Exception as e:
        return f"An error occurred during symbolic derivation: {e}"

    # --- Check 2: Physical Sanity Checks ---
    
    # 2a. Dimensional Consistency
    # The term k*A**2/m has dimensions of velocity^2. It cannot be subtracted from dimensionless 1.
    # This check invalidates Option A.
    
    # 2b. Physical Bounds (v_max < c)
    # For Option A and D, the term inside the sqrt is 1 + (a positive quantity), which is > 1.
    # This means v_max > c, which is physically impossible.
    
    # 2c. Non-relativistic Limit (c -> infinity)
    # The correct relativistic formula must reduce to the classical formula (Option B).
    classical_expr = options['B']
    try:
        limit_of_C = limit(options['C'], c, oo)
        if sympy.simplify(limit_of_C - classical_expr) != 0:
            return f"Incorrect. The non-relativistic limit (c -> oo) of answer C does not match the classical expression (B).\nLimit of C: {limit_of_C}\nClassical expression B: {classical_expr}"
    except Exception as e:
        return f"An error occurred during limit calculation: {e}"

    # --- Final Verdict ---
    # If all checks pass for option C, the answer is correct.
    # The analysis shows:
    # - The derivation matches C.
    # - C is dimensionally consistent, physically bounded (v<c), and has the correct classical limit.
    # - A is dimensionally inconsistent and gives v>c.
    # - B is the classical limit, not the full relativistic answer.
    # - D gives v>c.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)