import sympy

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio between the lightest solitonic
    excitation and the subsequent higher excitation in the CP(N-1) model.
    """
    # Define N as a symbolic variable for the limit calculation.
    N = sympy.Symbol('N')
    pi = sympy.pi

    # In the large-N CP(N-1) model, the mass M_k of a solitonic state with charge k
    # is proportional to sin(pi*k/N). The proportionality constant cancels out in the ratio.
    # The states are labeled by k = 1, 2, ..., N-1.

    # The lightest excitation corresponds to the smallest k, which is k=1.
    k_lightest = 1
    mass_lightest_term = sympy.sin(pi * k_lightest / N)

    # The subsequent higher excitation corresponds to the next value, k=2.
    k_next = 2
    mass_next_term = sympy.sin(pi * k_next / N)

    # The mass ratio is the ratio of the term for k=2 to the term for k=1.
    mass_ratio_expr = mass_next_term / mass_lightest_term

    # Calculate the asymptotic limit of this ratio as N approaches infinity.
    asymptotic_ratio = sympy.limit(mass_ratio_expr, N, sympy.oo)

    # Per the instructions, we print the final equation showing the numbers involved.
    # The equation shows the ratio of the mass for the state k=2 to the state k=1 equals the final result.
    print(f"Ratio of Mass(k={k_next}) to Mass(k={k_lightest}) = {asymptotic_ratio}")

if __name__ == '__main__':
    solve_mass_ratio()