import sympy

def solve_r0f_expression():
    """
    This function derives and prints the symbolic expression for R0f.
    """
    # Define the symbolic variables using names that match the problem description.
    # Greek letters are spelled out for variable names but represented with symbols for printing.
    b, pg, c, pt = sympy.symbols('b p_g c p_t')
    gamma_t = sympy.Symbol('𝛾_t')
    mu_t = sympy.Symbol('μ_t')
    mu_g = sympy.Symbol('μ_g')
    tau = sympy.Symbol('τ')

    # R0f is the total number of new tree fires caused by a single burning tree over its lifetime.
    # The infection path is Tree -> Grass -> Tree.

    # 1. Number of grass patches ignited by one tree over its lifetime.
    #    Rate of ignition = b * pg
    #    Lifetime of burning tree = 1 / (gamma_t + mu_t)
    #    Total grass ignited = (b * pg) / (gamma_t + mu_t)
    
    # 2. Number of trees ignited by one grass patch, considering its own lifecycle.
    #    Probability the grass survives latency to become infectious = exp(-mu_g * tau)
    #    Lifetime of an infectious grass patch = 1 / mu_g
    #    Rate of tree ignition from one infectious grass patch = c * pt
    #    Total trees ignited per infectious grass patch = (c * pt) / mu_g
    
    # R0f is the product of the number of grass fires created by a tree,
    # the probability of those grass fires becoming infectious, and the
    # number of new tree fires each infectious grass patch creates.
    
    # Combining the logic:
    # R0f = (Total grass ignited by 1 tree) * (Prob. grass becomes infectious) * (Total trees ignited per infectious grass)
    # This can be formulated as:
    # R0f = (Rate of creating infectious grass) * (Lifetime of tree) * (Number of trees per infectious grass)
    
    # Final derived expression:
    numerator = b * pg * c * pt * sympy.exp(-mu_g * tau)
    denominator = mu_g * (gamma_t + mu_t)
    
    R0f_expr = numerator / denominator

    # Print the final equation clearly, showing each component.
    print("The expression for R0f is:")
    # Using sympy.pretty provides a nicely formatted multiline output
    print(sympy.pretty(R0f_expr, use_unicode=True))
    
    # To satisfy the format requirement, we will also provide a single-line version
    # and the final answer in the specified format.
    final_expression_str = f"({b}*{pg}*{c}*{pt}*exp(-{mu_g}*{tau})) / ({mu_g}*({gamma_t} + {mu_t}))"

    # print("\nFinal expression in a single line:")
    # print("R0f =", final_expression_str)
    
solve_r0f_expression()