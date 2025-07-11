import sympy

def solve_yukawa_ratio():
    """
    Calculates the ratio R for the given Yukawa theory.

    The one-loop counter-terms in the MS-bar scheme are used.
    The calculation is based on the following key results from one-loop computations
    in D=4-2*epsilon dimensions, where the divergent part is proportional
    to C = g**2 / (16 * pi**2 * epsilon):

    1. From the fermion self-energy diagram:
       - The fermion wave function renormalization is delta_Zx = -C / 2.
       - The sum of fermion wave function and mass counter-terms is delta_Zx + delta_Zmx = C.

    2. From the vertex correction diagram (assuming delta_Z_phi = 0):
       - The sum of Yukawa coupling and fermion wave function counter-terms is
         delta_Zg + delta_Zx = -C.

    This script uses these relations to find delta_Zmx and delta_Zg in terms of C,
    and then computes the ratio R = delta_Zx / (delta_Zg + delta_Zmx).
    """

    # Define a symbolic constant C for the common factor g^2 / (16 * pi^2 * epsilon)
    # This factor will cancel out in the final ratio.
    C = sympy.Symbol('C')

    # From fermion self-energy: delta_Zx = -g^2 / (32 * pi^2 * epsilon)
    delta_Zx = -C / 2

    # From fermion self-energy: delta_Zx + delta_Zmx = g^2 / (16 * pi^2 * epsilon)
    # => delta_Zmx = C - delta_Zx
    delta_Zmx = C - delta_Zx

    # From vertex correction (given delta_Zphi = 0): delta_Zg + delta_Zx = -g^2 / (16 * pi^2 * epsilon)
    # => delta_Zg = -C - delta_Zx
    delta_Zg = -C - delta_Zx

    # Now, calculate the ratio R
    R = delta_Zx / (delta_Zg + delta_Zmx)

    # Let's define the numerical coefficients for clarity in the final output.
    # We express each counter-term as coeff * (g^2 / (32 * pi^2 * epsilon))
    unit_U = C / 2
    
    coeff_Zx = sympy.simplify(delta_Zx / unit_U)
    coeff_Zmx = sympy.simplify(delta_Zmx / unit_U)
    coeff_Zg = sympy.simplify(delta_Zg / unit_U)
    
    # Print the explanation and the final equation with numerical coefficients.
    print("The one-loop counter-term coefficients can be expressed in units of U = g^2 / (32 * pi^2 * epsilon):")
    print(f"delta_Zx = {coeff_Zx} * U")
    print(f"delta_Zmx = {coeff_Zmx} * U")
    print(f"delta_Zg = {coeff_Zg} * U")
    print("\nThe ratio R is defined as R = delta_Zx / (delta_Zg + delta_Zmx).")
    print("Substituting the coefficients:")
    
    # Print the equation with the numbers
    numerator_val = float(coeff_Zx)
    denominator_val = float(coeff_Zg) + float(coeff_Zmx)
    
    print(f"R = ({numerator_val}) / (({coeff_Zg}) + ({coeff_Zmx}))")
    print(f"R = {numerator_val} / {denominator_val}")
    
    # Print the final answer
    final_R = sympy.simplify(R)
    print(f"\nThe calculated value of the ratio R is: {final_R}")

solve_yukawa_ratio()
<<< -1/2 >>>