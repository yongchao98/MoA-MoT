import sympy
from sympy import sin, pi, limit, cos

def solve_mass_ratio():
    """
    Calculates the asymptotic mass ratio for excitations in the CP(N-1) model.
    """
    # Define symbols for symbolic computation
    N = sympy.Symbol('N', integer=True, positive=True)

    # Step 1: Explain the theoretical background and the mass formula.
    print("In the two-dimensional CP(N-1) sigma model, the spectrum of stable particles is known from exact S-matrix theory.")
    print("The masses of these particles, which can be considered solitonic excitations, follow a specific pattern.")
    print("The mass of the k-th excitation (m_k) relative to the lightest one (m_1) is given by:")
    print("m_k / m_1 = sin(k*pi/N) / sin(pi/N)\n")

    # Step 2: Identify the particles for the ratio.
    # The lightest excitation is the k=1 particle. Its mass is m_1.
    # The subsequent higher excitation is the k=2 particle.
    print("The 'lightest solitonic excitation' corresponds to k=1.")
    print("The 'subsequent higher excitation' corresponds to k=2.\n")

    # Step 3: Formulate the mass ratio R.
    print("The mass ratio R is therefore m_2 / m_1.")
    # The expression for the ratio R(N)
    ratio_expr = sin(2 * pi / N) / sin(pi / N)
    print(f"R(N) = sin(2*pi/N) / sin(pi/N)")

    # Simplify the expression using trigonometric identities
    # sin(2x) = 2*sin(x)*cos(x)
    ratio_expr_simplified = 2 * cos(pi / N)
    print(f"Using sin(2x) = 2*sin(x)*cos(x), this simplifies to: R(N) = {ratio_expr_simplified}\n")

    # Step 4: Calculate the asymptotic limit as N approaches infinity.
    print("To find the asymptotic ratio, we compute the limit of R(N) as N -> infinity.")
    asymptotic_ratio = limit(ratio_expr_simplified, N, sympy.oo)

    # Step 5: Print the final calculation and result.
    print("The final equation is the evaluation of this limit:")
    # The numbers in the final equation are 2 and the result of cos(0) which is 1.
    numerator = 2
    denominator = 1
    final_result = asymptotic_ratio
    print(f"lim_{{N->oo}} (2 * cos(pi/N)) = {numerator} * cos(0) = {numerator} * {denominator} = {final_result}")
    print(f"\nThe leading-order asymptotic mass ratio is: {final_result}")

if __name__ == '__main__':
    solve_mass_ratio()