import sympy

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio between the lightest
    and subsequent higher solitonic excitations in the CP(N-1) model.
    """

    # 1. Define the quantum numbers for the states in question.
    # The "lightest solitonic excitation" corresponds to the fundamental kink state, k=1.
    k_lightest = 1
    # The "subsequent higher excitation" corresponds to the first bound state of kinks, k=2.
    k_higher = 2

    # 2. Formulate the mass ratio expression.
    # The mass of a kink bound state 'k' is proportional to sin(k*pi/N).
    # The proportionality constant cancels out in the ratio.
    N = sympy.Symbol('N')
    mass_ratio_expr = sympy.sin(k_higher * sympy.pi / N) / sympy.sin(k_lightest * sympy.pi / N)

    # 3. Compute the limit as N approaches infinity.
    # As N -> infinity, let x = pi/N, then x -> 0.
    # The expression becomes sin(2x) / sin(x) = (2*sin(x)*cos(x)) / sin(x) = 2*cos(x).
    # The limit is 2*cos(0) = 2.
    asymptotic_ratio = sympy.limit(mass_ratio_expr, N, sympy.oo)

    # 4. Output the final result and the numbers in the final equation.
    print("The mass ratio is given by the equation: M_k2 / M_k1 = R")
    print("Where k1 is the quantum number of the lightest excitation and k2 is that of the next.")
    print("\nThe final equation with its numerical values is:")
    print(f"Mass(k={k_higher}) / Mass(k={k_lightest}) = {int(asymptotic_ratio)}")

    print("\nThe individual numbers in this final equation are:")
    print(k_higher)
    print(k_lightest)
    print(int(asymptotic_ratio))

solve_mass_ratio()