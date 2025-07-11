import math

def solve():
    """
    This function derives the analytical approximation for the integral I(epsilon)
    for small epsilon.
    """
    # Parameters from the dominant term in the denominator polynomial (9.0*x^5.0).
    c_min = 9.0
    p_min = 5.0

    # The approximation has the form I(epsilon) ~ C * epsilon^a.
    
    # Calculate the exponent 'a'.
    a = (1.0 / p_min) - 1.0

    # Calculate the coefficient 'C'.
    # The formula for C is: pi / (p_min * c_min^(1/p_min) * sin(pi/p_min)).
    pi = math.pi
    c_min_term = c_min**(1.0 / p_min)
    sin_term = math.sin(pi / p_min)
    C = pi / (p_min * c_min_term * sin_term)

    # Print the step-by-step derivation of the formula with numerical values.
    print("The analytical formula that approximates I(epsilon) for small epsilon is of the form: C * epsilon^a.")
    print("\nDerivation based on approximation near x=0:")
    print(f"1. Approximate the denominator term f(x) by its lowest power term: {c_min} * x^{p_min}")
    print("2. Approximate the integral as I(e) ≈ Integral[0, inf] dx / (e + c_min * x^p_min).")
    print("3. This integral evaluates to C * e^a, where:")
    print(f"   a = (1 / p_min) - 1")
    print(f"   C = pi / (p_min * c_min^(1/p_min) * sin(pi/p_min))")

    print("\nPlugging in the values from the problem (c_min=9.0, p_min=5.0):")
    print(f"a = (1 / {p_min}) - 1 = {a:.1f}")
    
    # Show the formula for C with numbers.
    formula_c_str = f"C = {pi:.6f} / ({p_min} * {c_min}^(1/{p_min}) * sin({pi:.6f}/{p_min}))"
    print(formula_c_str)
    print(f"C = {pi:.6f} / ({p_min} * {c_min_term:.6f} * {sin_term:.6f})")
    print(f"C = {pi:.6f} / {(p_min * c_min_term * sin_term):.6f}")
    print(f"C ≈ {C:.6f}")

    print("\n------------------------------------------------------")
    print("Final Analytical Formula:")
    # Print the formula with calculated numeric values for C and a.
    print(f"I(epsilon) ≈ {C:.4f} * epsilon^({a:.1f})")
    
    print("\nExplicit formula showing all the numerical constants:")
    # Print the final equation with all numbers in it, as requested.
    print(f"I(epsilon) ≈ ({pi:.4f} / ({p_min} * {c_min}**(1.0/{p_min}) * math.sin({pi:.4f}/{p_min}))) * epsilon**((1.0/{p_min}) - 1.0)")
    print("------------------------------------------------------")

solve()