import sympy
from sympy import Symbol, Rational, integrate, log, exp, latex

def solve_queueing_problem():
    """
    Solves the queueing theory problem described by the user.
    1. Sets up the expression for the mean response time x.
    2. Symbolically computes the integral to find x.
    3. Filters the terms of x based on the given rules.
    4. Prints the result in LaTeX format.
    """
    # Define symbols and constants
    s = Symbol('s')
    t = Symbol('t')
    lambda_val = Rational(3, 2)

    # The mean response time for a job of size s is E[T(s)].
    # E[T(s)] = ( (lambda * Integral(t**2, 0, s)) / 2 + s ) / (1 - lambda * Integral(t, 0, s))
    
    # Let's calculate the components
    # rho_s = lambda * Integral(t dt) from 0 to s
    rho_s = lambda_val * integrate(t, (t, 0, s))
    
    # Numerator of E[T(s)]
    numerator_E_Ts = (lambda_val * integrate(t**2, (t, 0, s))) / 2 + s
    
    # Denominator of E[T(s)]
    denominator_E_Ts = 1 - rho_s
    
    # Full expression for E[T(s)]
    E_Ts = numerator_E_Ts / denominator_E_Ts
    E_Ts_simplified = sympy.simplify(E_Ts)
    
    # The overall mean response time, x, is the integral of E[T(s)] from 0 to 1
    x = integrate(E_Ts_simplified, (s, 0, 1))
    
    print(f"The formula for the conditional mean response time E[T(s)] simplifies to: {E_Ts_simplified}")
    print(f"Integrating from s=0 to s=1, the optimal mean response time x is: {x}\n")

    # The problem asks to remove certain terms from x.
    # Let's decompose x into its additive terms.
    if isinstance(x, sympy.Add):
        terms = x.args
    else:
        terms = [x]

    print("Analyzing the terms of x:")
    
    remaining_terms = []
    
    for term in terms:
        is_rational_term = False
        is_log_of_rational_term = False
        
        # Check 1: Is the term a rational number?
        if term.is_rational:
            is_rational_term = True
            print(f"- Term '{term}' is a rational number. It will be removed.")
        
        # Check 2: Is the term a logarithm of a rational number?
        # This is true if exp(term) is a rational number.
        else:
            exp_term = exp(term)
            # Use sympy.ask to check if the result is rational, handling symbolic expressions
            if sympy.ask(sympy.Q.rational(exp_term)):
                is_log_of_rational_term = True
                print(f"- Term '{term}': exp({term}) is {exp_term}, which is rational. It will be removed.")
            else:
                # Add the simplified version of the exponentiation for clarity
                exp_term_simplified = sympy.simplify(exp_term)
                print(f"- Term '{term}': exp({term}) is {exp_term_simplified}, which is not rational. It remains.")

        if not is_rational_term and not is_log_of_rational_term:
            remaining_terms.append(term)
            
    # The final answer is the sum of the remaining terms.
    final_answer_expr = sum(remaining_terms)
    
    # Format the final answer in LaTeX
    # We build the string manually for better formatting, e.g., \frac{a}{b} instead of a/b
    if final_answer_expr.is_Mul:
        c, m = final_answer_expr.as_coeff_Mul()
        if m.func == log:
            arg = m.args[0]
            if c.is_Rational:
                latex_str = f"\\frac{{{c.p}}}{{{c.q}}} \\ln{{{latex(arg)}}}"
            else:
                 latex_str = latex(final_answer_expr) # fallback
        else:
            latex_str = latex(final_answer_expr) # fallback
    else:
        latex_str = latex(final_answer_expr)
    
    print(f"\nThe remaining term of x is: {final_answer_expr}")
    print(f"This formatted in standard LaTeX notation is: {latex_str}")

if __name__ == '__main__':
    solve_queueing_problem()
