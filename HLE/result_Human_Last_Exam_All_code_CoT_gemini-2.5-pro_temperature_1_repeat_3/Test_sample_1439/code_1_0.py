import sympy as sp

def solve_critical_exponent_order():
    """
    Symbolically derives the order of the first correction to the critical exponent nu
    in phi^4 theory using the perturbative renormalization group framework.
    """
    # 1. Define the symbolic variables
    # u: the renormalized coupling constant
    # epsilon: the deviation from the upper critical dimension (d=4), i.e., epsilon = 4 - d
    # A, B: positive constants from loop calculations, their specific values are not needed
    #         to determine the order of the correction.
    u, epsilon, A, B = sp.symbols('u epsilon A B', positive=True)

    # 2. Define the fundamental RG equations to the lowest non-trivial order (one-loop).
    # The beta function describes the scale dependence of the coupling u.
    beta_u = -epsilon * u + A * u**2

    # The anomalous dimension of the phi^2 operator, gamma_2(u), is related to nu.
    # Its first contribution comes from a one-loop diagram and is linear in u.
    gamma2_u = B * u

    # 3. The critical exponent nu is related to gamma_2 via nu^{-1} = 2 - gamma_2(u).
    # We can express nu as a function of the coupling u.
    nu_of_u = 1 / (2 - gamma2_u)

    # 4. The question asks for the order in 'u' of the first correction.
    # To find this, we can perform a Taylor series expansion of nu(u) around u=0.
    # We expand to the second order to see the constant term and the first correction.
    nu_series_in_u = nu_of_u.series(u, 0, 2)

    # 5. Print the derivation and the final result clearly.
    print("This script determines the order of the first correction to the critical exponent nu.")
    print("The derivation follows the standard perturbative epsilon-expansion for phi^4 theory.")
    print("-" * 75)

    print("Step 1: Define the core RG equations.")
    print(f"The beta function is: beta(u) = {beta_u}")
    print(f"The anomalous dimension is: gamma_2(u) = {gamma2_u}")
    print("-" * 75)

    print("Step 2: Express the critical exponent nu in terms of the coupling u.")
    print(f"The relation is: nu(u) = 1 / (2 - gamma_2(u))")
    print(f"Substituting gamma_2 gives: nu(u) = {nu_of_u}")
    print("-" * 75)

    print("Step 3: Perform a series expansion of nu(u) around u=0 to find the corrections.")
    print("The mean-field value of nu is the constant term (at u=0).")
    print("The first correction is the next term in the series.")
    print("\nThe series expansion is:")
    print(f"nu(u) = {nu_series_in_u}")
    print("-" * 75)

    # Extract the order and the coefficient of the first correction term.
    # The coefficient of u in the series expansion.
    correction_coeff = nu_series_in_u.coeff(u, 1)
    # The order of u is 1.
    correction_order = 1

    print("Final Analysis:")
    print(f"The mean-field value (the u^0 term) is 1/2.")
    print(f"The first non-vanishing contribution (correction) to this value is the term '{correction_coeff}*u'.")
    print("\nThis contribution is of order 1 in the coupling constant u.")
    print("\nFinal equation for nu with the first correction:")
    # The "final equation" showing the numbers as requested by the prompt.
    final_eq_str = f"nu = 1/2 + ({B}/4)*u + O(u^2)"
    print(final_eq_str)
    print("The numbers in this equation are: 1, 2, 4, 2.")

if __name__ == '__main__':
    solve_critical_exponent_order()

<<<1>>>