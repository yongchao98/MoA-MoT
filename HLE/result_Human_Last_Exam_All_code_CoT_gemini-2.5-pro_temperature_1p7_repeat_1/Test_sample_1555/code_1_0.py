import sympy

def calculate_cpn1_mass_ratio():
    """
    Determines the leading-order asymptotic mass ratio between the lightest (k=1)
    and subsequent higher (k=2) solitonic excitations in the CP(N-1) model.

    This is achieved by taking the limit of the mass ratio M_2 / M_1 as N -> infinity,
    using the semi-classical mass formula M_k ∝ sin(pi*k/N).
    """
    # Define N as a symbolic variable for the limit calculation
    N = sympy.Symbol('N')

    # The excitations are labeled by an integer k.
    # The lightest solitonic excitation corresponds to k=1.
    k1 = 1
    # The subsequent higher excitation corresponds to k=2.
    k2 = 2

    print(f"The physical model is the CP(N-1) sigma model in the large-N limit.")
    print(f"The mass of the k-th excitation is proportional to sin(pi*k/N).")
    print("-" * 50)

    # Proportional expression for the mass of the lightest excitation (k=1)
    # The constant of proportionality cancels in the ratio.
    m1_expr = sympy.sin(sympy.pi * k1 / N)

    # Proportional expression for the mass of the next excitation (k=2)
    m2_expr = sympy.sin(sympy.pi * k2 / N)

    # The ratio of the masses, M_2 / M_1
    ratio_expr = m2_expr / m1_expr

    print(f"We are calculating the ratio between the mass of the k={k2} state and the k={k1} state.")
    print(f"The equation for the mass ratio is:")
    print(f"Ratio = M_{k2} / M_{k1} = sin({k2}*pi/N) / sin({k1}*pi/N)")
    print("-" * 50)
    
    # Calculate the asymptotic limit as N approaches infinity
    asymptotic_ratio = sympy.limit(ratio_expr, N, sympy.oo)

    print(f"To find the asymptotic value, we compute the limit as N approaches infinity:")
    print(f"lim (N->inf) [ {ratio_expr} ] = {asymptotic_ratio}")
    print("-" * 50)
    print(f"The final calculated number for the asymptotic mass ratio is {asymptotic_ratio}.")

if __name__ == '__main__':
    calculate_cpn1_mass_ratio()

<<<2>>>