import math

def solve_integral_approximation():
    """
    Calculates and prints the analytical approximation for the integral I(epsilon)
    for the small epsilon regime.
    """
    # From the analysis, the integral I(epsilon) is approximated by C * epsilon^p.
    
    # The exponent 'p' is derived from the dominant term's power 'n'.
    # g(x) is dominated by 9x^5, so n=5.
    n = 5.0
    p = (1.0 / n) - 1.0
    
    # The coefficient 'C' depends on the coefficient of the dominant term 'a' and 'n'.
    # The dominant term is a*x^n = 9*x^5.
    a = 9.0
    
    # The coefficient C = a^(-1/n) * integral from 0 to inf of 1/(1+u^n) du
    # The integral is known to be (pi/n) / sin(pi/n)
    integral_value = (math.pi / n) / math.sin(math.pi / n)
    C = (a**(-1.0/n)) * integral_value
    
    # Output the results, including each number in the final equation.
    print("The analytical approximation for the integral is of the form: I(epsilon) = C * epsilon^p")
    print("\n--- Derivation of the components ---")
    
    # Exponent p
    print(f"The exponent p is calculated as 1/n - 1, where n is the lowest power of x (n=5.0).")
    print(f"p = 1.0/{n} - 1.0 = {p}")
    
    # Coefficient C
    print("\nThe coefficient C is calculated as a^(-1/n) * [pi / (n * sin(pi/n))], where a=9.0 and n=5.0.")
    a_inv_n = a**(-1.0/n)
    pi_over_n_sin = (math.pi / n) / math.sin(math.pi / n)
    print(f"a^(-1/n) = {a}^(-1.0/{n}) = {a_inv_n}")
    print(f"[pi / (n * sin(pi/n))] = [pi / ({n} * sin(pi/{n}))] = {pi_over_n_sin}")
    print(f"C = {a_inv_n} * {pi_over_n_sin} = {C}")
    
    # Final formula
    print("\n--- Final Approximate Formula ---")
    print(f"I(epsilon) approx {C} * epsilon^({p})")

solve_integral_approximation()