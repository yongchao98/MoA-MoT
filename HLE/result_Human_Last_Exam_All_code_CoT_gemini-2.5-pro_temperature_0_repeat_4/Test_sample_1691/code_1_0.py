import math

def solve_integral_approximation():
    """
    This function derives and prints the analytical approximation for the given integral
    in the small epsilon regime.
    """
    
    # The integral is I(epsilon) = integral from 0 to 15 of 1 / (epsilon + g(x)) dx
    # where g(x) = 9.0 * x**5.0 + 5.0 * x**6.0 + 9.0 * x**8.0

    # For small epsilon, the integral is dominated by the behavior of g(x) near x=0.
    # We approximate g(x) by its leading term, g(x) ~ a * x**n.
    a = 9.0
    n = 5.0

    # The approximation for the integral takes the form: I(epsilon) ~ C * epsilon**p
    
    # The power p is determined by the scaling and is given by (1-n)/n.
    p = (1.0 - n) / n

    # The coefficient C is given by a**(-1/n) * integral_0^inf(1/(1+u**n)du)
    # The integral part is equal to (pi/n) / sin(pi/n).
    
    # Let's calculate the components of C.
    # First component from the coefficient 'a'
    term1 = a**(-1.0/n)
    
    # Second component from the definite integral
    term2 = (math.pi / n) / math.sin(math.pi / n)
    
    # The final coefficient C
    C = term1 * term2

    # Now, we print the derivation and the final formula.
    print("The analytical approximation for the integral I(epsilon) for small epsilon is derived as follows:")
    print("1. The integral is dominated by the region near x=0.")
    print(f"2. The function in the denominator is approximated by its leading term: g(x) ≈ {a} * x^{n}")
    print("3. The integral is approximated by I(epsilon) ≈ integral from 0 to infinity of 1 / (epsilon + a * x^n) dx.")
    print("4. This leads to an asymptotic formula of the form: I(epsilon) ≈ C * epsilon^p")
    
    print("\nCalculating the values for p and C:")
    print(f"The power p = (1 - n) / n = (1 - {n}) / {n} = {p}")
    print(f"The coefficient C = a^(-1/n) * (pi/n) / sin(pi/n)")
    print(f"C = ({a})^(-1/{n}) * (pi/{n}) / sin(pi/{n})")
    print(f"C = {term1:.6f} * {term2:.6f}")
    print(f"C = {C:.6f}")

    print("\nThus, the final analytical formula is:")
    print(f"I(epsilon) ≈ {C:.4f} * epsilon^({p})")

solve_integral_approximation()