import sympy as sp

def solve_invariant_density():
    """
    Calculates the normalised density of the invariant measure for the map
    T(x) = 1/sqrt(x) mod 1.

    The functional form of the density is known to be rho(x) = C*(1 + 2*sqrt(x)).
    This script computes the normalization constant C and prints the final result.
    """

    # Define the symbolic variable and constant
    x = sp.Symbol('x', real=True, positive=True)
    C = sp.Symbol('C')

    # Define the unnormalized density function
    # Using x**(1/2) for sqrt(x)
    rho_unnormalized = C * (1 + 2 * x**sp.S.Half)

    # Calculate the integral of the unnormalized density from 0 to 1
    integral_val = sp.integrate(rho_unnormalized, (x, 0, 1))

    # For a normalized density, the integral must equal 1.
    # Solve for the normalization constant C.
    normalization_equation = sp.Eq(integral_val, 1)
    C_val = sp.solve(normalization_equation, C)[0]

    # Construct the final normalized density function
    rho_normalized = rho_unnormalized.subs(C, C_val)
    
    # Print the results as requested
    print("The functional form of the invariant density is rho(x) = C * (1 + 2*sqrt(x)).")
    print("To find the normalization constant C, we solve the equation: Integral(rho(x) dx) from 0 to 1 = 1.")
    print(f"The symbolic integral evaluates to: {integral_val.subs(C, 'C')}")
    print(f"Solving '{integral_val.subs(C, 'C')} = 1' gives C = {C_val}.")
    print("\n--------------------------------------------------")
    print("The final expression for the normalised density is:")
    
    # To satisfy the output format of showing each number
    num, den = C_val.as_numer_denom()
    expr_str = f"rho(x) = ({num}/{den}) * (1 + 2 * x**0.5)"
    print(expr_str)

    print("\nThe numbers in the final equation are:")
    print(f"  - The numerator of the leading constant: {num}")
    print(f"  - The denominator of the leading constant: {den}")
    print(f"  - The additive term inside the parenthesis: 1")
    print(f"  - The coefficient of the sqrt term inside the parenthesis: 2")
    print(f"  - The power of x in the sqrt term: 0.5")
    print("--------------------------------------------------")


if __name__ == '__main__':
    solve_invariant_density()
