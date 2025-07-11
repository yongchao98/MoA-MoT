import sympy as sp

def solve_limit():
    """
    Symbolically calculates the limit based on the derived asymptotic behavior of p_n.
    """
    n = sp.Symbol('n', real=True, positive=True)
    # C is an unknown positive constant from the asymptotic relation p_n ~ C/n
    C = sp.Symbol('C', real=True, positive=True)

    # We derived that p_n is proportional to 1/n for large n
    p_n = C / n
    
    # The expression to be evaluated
    expression = sp.log(1 / p_n) / sp.log(n)
    
    # Calculate the limit as n tends to infinity
    limit_result = sp.limit(expression, n, sp.oo)
    
    # Print the steps of the derivation as requested
    print("Let the escape probability p_n have the asymptotic behavior p_n = C/n for some constant C.")
    print("The quantity to find is the limit of ln(1/p_n) / ln(n) as n -> infinity.\n")
    print("Step 1: Express 1/p_n")
    print(f"1/p_n = 1 / ({p_n}) = {1/p_n}\n")
    
    print("Step 2: Take the natural logarithm")
    # Use expand_log to show the separation
    ln_1_over_pn = sp.expand_log(sp.log(1/p_n), force=True)
    print(f"ln(1/p_n) = ln({1/p_n}) = {ln_1_over_pn}\n")

    print("Step 3: Form the full expression")
    full_expression_latex = f"({sp.latex(ln_1_over_pn)}) / ({sp.latex(sp.log(n))})"
    print(f"The expression is: {full_expression_latex}")
    simplified_expression = sp.simplify(expression)
    print(f"This simplifies to: {simplified_expression}\n")

    print("Step 4: Calculate the limit")
    print(f"lim_{{n->oo}} ({simplified_expression}) = 1 - 0 = 1\n")

    print("Final equation with numbers:")
    print(f"lim (1 - ln(C)/ln(n)) = 1 - 0 = 1")

solve_limit()
