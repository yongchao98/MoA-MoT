import sympy

def print_overlap_integral_2s():
    """
    This function prints the derived analytical expression for the 2s-2s overlap
    integral in the H2+ ion.
    """
    # Define the symbolic variables
    R = sympy.Symbol('R', real=True, positive=True)
    zeta = sympy.Symbol('ζ', real=True, positive=True) # Using 'zeta' for the Greek letter ζ
    K = sympy.Symbol('K')

    # Define the polynomial part of the solution
    # S(K) = exp(-K) * [1 + K + (4/5)K^2 + (1/3)K^3 + (1/15)K^4]
    polynomial_part = 1 + K + (sympy.Rational(4, 5) * K**2) + \
                      (sympy.Rational(1, 3) * K**3) + (sympy.Rational(1, 15) * K**4)

    # Construct the full expression for the overlap integral S
    S = sympy.exp(-K) * polynomial_part

    # Print the final formula in a readable format
    print("The overlap integral S for two 2s orbitals in H₂⁺ is given by:")
    print("-" * 60)
    print(f"S(R, ζ) = {S}\n")
    print("Where K is defined as:")
    print(f"K = ({zeta} * R) / 2")
    print("-" * 60)

    # The problem asks to output each number in the final equation.
    # We will manually construct the string to match the requested format.
    # S = e^(-K) * (1 + 1*K + 4/5*K^2 + 1/3*K^3 + 1/15*K^4)
    eq_str = (
        "S = exp(-K) * (1 + K + (4/5)*K**2 + (1/3)*K**3 + (1/15)*K**4)"
    )
    print("\nFinal equation with explicit coefficients:")
    print(eq_str)


if __name__ == "__main__":
    print_overlap_integral_2s()