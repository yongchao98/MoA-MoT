import sympy

def polymer_force_law():
    """
    Calculates and prints the symbolic expression for the force law of a
    thermally isolated polymer chain.
    """
    # Define the symbolic variables for the equation
    # x: separation of the polymer ends
    # ell: length of each strut (monomer)
    # n: number of segments
    # E0: kinetic energy of the polymer at zero extension (x=0)
    x = sympy.Symbol('x')
    ell = sympy.Symbol('ell')
    n = sympy.Symbol('n', integer=True, positive=True)
    E0 = sympy.Symbol('E(0)', positive=True)

    # From the derivation, the energy E as a function of extension x is:
    # E(x) = E(0) * exp(x^2 / (n^2 * l^2))
    # We use the large n approximation (3n-1 -> 3n).
    energy_expr = E0 * sympy.exp(x**2 / (n**2 * ell**2))

    # The force F is the derivative of the energy E with respect to x.
    force_expr = sympy.diff(energy_expr, x)

    # The resulting force law has the form:
    # F(x) = A * x * exp(B * x**2)
    # We will identify and print the components A and B.

    # Extract the pre-exponential factor (coefficient of x)
    prefactor = 2 * E0 / (n**2 * ell**2)
    
    # Extract the coefficient of x^2 in the exponent
    exp_coefficient = 1 / (n**2 * ell**2)

    print("The force law F(x) for a thermally isolated polymer chain under small extension x is:")
    print(f"F(x) = {force_expr}\n")

    print("This can be written in the form: F(x) = (A) * x * exp(B * x**2)")
    print("where the coefficients are:\n")
    
    print("Coefficient A (Prefactor):")
    # To satisfy the "output each number" requirement, let's break it down
    print(f"  Numerator: 2 * E(0)")
    print(f"  Denominator: n**2 * ell**2")
    print(f"  A = {prefactor}\n")

    print("Coefficient B (in the exponent):")
    print(f"  Numerator: 1")
    print(f"  Denominator: n**2 * ell**2")
    print(f"  B = {exp_coefficient}\n")


if __name__ == '__main__':
    polymer_force_law()