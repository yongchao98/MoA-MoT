import sympy

def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4) without X
    for degrees up to 100.

    The Poincare series for the equivariant cohomology A = H_{SO(4)}^*(SO(4) \setminus X) is given by
    P_A(t) = (1 + t^3 + 2*t^6) / (1 - t^4)^2.
    We need to find the sum of the coefficients of the Taylor expansion of this series
    for powers of t from 0 to 100.
    """
    t = sympy.Symbol('t')
    
    # Numerator of the Poincare series
    numerator = 1 + t**3 + 2*t**6
    
    # Denominator of the Poincare series
    denominator = (1 - t**4)**2
    
    # The Poincare series as a rational function
    poincare_series = numerator / denominator
    
    # We want to find the coefficients of the series expansion up to degree 100.
    # The 'series' method in sympy expands the function around a point (here t=0).
    # n=101 to include the coefficient of t^100.
    series_expansion = poincare_series.series(t, 0, 101).removeO()
    
    # Get the dictionary of coefficients. The keys are the powers of t, values are the coefficients.
    coeffs = series_expansion.as_poly().all_coeffs()
    
    # The coefficients are returned in descending order of power. We need to reverse them
    # to match ascending powers of t.
    coeffs.reverse()

    # We need to sum the coefficients a_k for k from 0 to 100.
    # We have to be careful with the list of coefficients. If the highest power in the
    # expansion is less than 100, we should only sum the available ones.
    # sympy's `all_coeffs` might not include zero coefficients for missing terms.
    # It is better to use `coeff` method for each power.
    
    total_rank = 0
    # Let's print out the calculation as in the plan for verification
    # a_{4j} = j+1
    # a_{4j+1} = 0
    # a_{4j+2} = 2j (for j>=1, a_2=0)
    # a_{4j+3} = j+1

    # sum for a_{4j}
    # 4j <= 100 => j <= 25
    sum_a_4j = 0
    for j in range(26):
        term_rank = j + 1
        print(f"rank of degree {4*j}: {term_rank}")
        sum_a_4j += term_rank
    
    # sum for a_{4j+1}, always 0
    
    # sum for a_{4j+2}
    # 4j+2 <= 100 => 4j <= 98 => j <= 24.5 => j <= 24
    sum_a_4j_2 = 0
    for j in range(25):
        if 4 * j + 2 > 100:
            continue
        term_rank = 0
        if j >= 1:
            term_rank = 2*j
        print(f"rank of degree {4*j + 2}: {term_rank}")
        sum_a_4j_2 += term_rank

    # sum for a_{4j+3}
    # 4j+3 <= 100 => 4j <= 97 => j <= 24.25 => j <= 24
    sum_a_4j_3 = 0
    for j in range(25):
        term_rank = j + 1
        print(f"rank of degree {4*j + 3}: {term_rank}")
        sum_a_4j_3 += term_rank
        
    total_rank = sum_a_4j + sum_a_4j_2 + sum_a_4j_3

    print(f"\nTotal rank for degrees <= 100:")
    print(total_rank)

solve()