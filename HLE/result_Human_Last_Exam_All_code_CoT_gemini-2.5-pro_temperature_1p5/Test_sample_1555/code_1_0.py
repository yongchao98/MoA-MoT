import sympy

def solve_mass_ratio():
    """
    Calculates the leading-order asymptotic mass ratio for excitations
    in the CP(N-1) model in the large-N limit.
    """
    # Step 1: Define the symbolic variables
    N = sympy.Symbol('N', positive=True, integer=True)
    m_g = sympy.Symbol('m_g', positive=True, real=True, nonzero=True) # Mass gap
    pi = sympy.pi

    # Step 2: Define the mass formula for the k-th excitation
    # M_k = (N * m_g / pi) * sin(k * pi / N)

    # The lightest solitonic excitation corresponds to k=1
    M1_expr = (N * m_g / pi) * sympy.sin(1 * pi / N)

    # The subsequent higher excitation corresponds to k=2
    M2_expr = (N * m_g / pi) * sympy.sin(2 * pi / N)

    # Step 3: Calculate the asymptotic mass for each excitation by taking the limit as N -> infinity
    # For k=1
    m1_asymptotic = sympy.limit(M1_expr, N, sympy.oo)
    # For k=2
    m2_asymptotic = sympy.limit(M2_expr, N, sympy.oo)
    
    # Step 4: The final ratio is m2_asymptotic / m1_asymptotic
    # Sympy will compute this as (2 * m_g) / (1 * m_g) = 2
    final_ratio = m2_asymptotic / m1_asymptotic
    
    # Extract the numerical coefficients for the final equation output
    # m1_asymptotic is of the form C1 * m_g, m2_asymptotic is C2 * m_g
    # We want to display the equation C2 / C1 = Ratio
    c1 = sympy.collect(m1_asymptotic, m_g).coeff(m_g)
    c2 = sympy.collect(m2_asymptotic, m_g).coeff(m_g)

    # Step 5: Print the results as per the requested format.
    # The final equation is Ratio = M2_limit / M1_limit
    print("In the large-N limit:")
    print(f"The asymptotic mass of the lightest excitation (k=1) is: M_1 = {c1} * m_g")
    print(f"The asymptotic mass of the next excitation (k=2) is: M_2 = {c2} * m_g")
    print("\nThe equation for the leading-order asymptotic mass ratio (M_2 / M_1) is:")
    # The prompt requests to output each number in the final equation.
    print(f"{int(c2)} / {int(c1)} = {int(final_ratio)}")

if __name__ == '__main__':
    solve_mass_ratio()