import sympy

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio between the lightest
    and subsequent higher excitation in the CP(N-1) model.
    """

    # Define N as a symbolic variable representing the number of colors,
    # approaching infinity.
    N = sympy.Symbol('N', real=True, positive=True)

    # In the large-N CP(N-1) model, the mass of the k-th bound state (solitonic
    # excitation) M_k is given by M_k = C * sin(k*pi/N), where C is a constant.
    # The lightest excitation corresponds to k=1.
    # The subsequent higher excitation corresponds to k=2.

    # We can write the expressions for their masses, omitting the constant C
    # as it will cancel out in the ratio.
    M1_expr = sympy.sin(sympy.pi / N)
    M2_expr = sympy.sin(2 * sympy.pi / N)

    # The mass ratio is M_2 / M_1.
    mass_ratio_expr = M2_expr / M1_expr

    # The problem asks for the asymptotic ratio as N -> infinity. We compute this
    # using a limit.
    limit_value = sympy.limit(mass_ratio_expr, N, sympy.oo)

    # To fulfill the request of showing the numbers in the final equation,
    # we can demonstrate the simplification and evaluation step-by-step.
    # The derivation uses the identity sin(2x) = 2*sin(x)*cos(x).
    # Ratio = sin(2*pi/N) / sin(pi/N)
    #       = (2 * sin(pi/N) * cos(pi/N)) / sin(pi/N)
    #       = 2 * cos(pi/N)
    # As N -> infinity, pi/N -> 0.
    # The limit becomes 2 * cos(0), which is 2 * 1.

    print("Step 1: Express the mass ratio.")
    print("The mass M_k is proportional to sin(k*π/N).")
    print("The ratio between the second (k=2) and first (k=1) excitation is:")
    print("Ratio = M_2 / M_1 = sin(2*π/N) / sin(π/N)\n")

    print("Step 2: Simplify the ratio using the double-angle identity sin(2x) = 2*sin(x)*cos(x).")
    print("Ratio = (2*sin(π/N)*cos(π/N)) / sin(π/N) = 2*cos(π/N)\n")

    print("Step 3: Calculate the asymptotic value by taking the limit as N approaches infinity.")
    print("lim (N->∞) [2*cos(π/N)] = 2 * cos(0)")
    
    # Final equation with numbers
    final_equation = f"2 * 1 = {int(limit_value)}"
    print(final_equation)


if __name__ == "__main__":
    solve_mass_ratio()