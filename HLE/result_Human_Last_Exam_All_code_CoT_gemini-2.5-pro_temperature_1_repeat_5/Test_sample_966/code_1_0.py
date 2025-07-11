import math

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q)
    for a complete intersection X of two quadrics in CP^102.
    """
    # Parameters from the problem statement
    n_ambient = 102
    degrees = [2, 2]

    # Step 1: Determine the dimension of the variety X
    dim_X = n_ambient - len(degrees)

    # Step 2: Determine the degree of X
    deg_X = 1
    for d in degrees:
        deg_X *= d

    # Step 3: Calculate the coefficient A = [h^m] in the formula for chi(X)
    # The term is (1+h)^(n+1) / product(1+d*h), where n=n_ambient, m=dim_X
    # A = [h^100] (1+h)^103 / (1+2h)^2
    # This is calculated by summing the convolution of the two series' coefficients:
    # A = sum_{j=0 to m} ( [h^{m-j}](1+h)^(n+1) ) * ( [h^j](1+2h)^-2 )
    # where [h^k](1+2h)^-2 = comb(-2, k) * 2^k = (-1)^k * (k+1) * 2^k
    
    N = n_ambient + 1
    m_coeff = dim_X
    A = 0
    
    # We are calculating A = sum_{j=0}^{100} C(103, 100-j) * (-1)^j * (j+1) * 2^j
    for j in range(m_coeff + 1):
        # Coefficient from (1+h)^103
        term1 = math.comb(N, m_coeff - j)
        
        # Coefficient from (1+2h)^-2
        term2 = ((-1)**j) * (j + 1) * (2**j)
        
        A += term1 * term2

    # Step 4: Calculate the Euler characteristic chi(X) = deg(X) * A
    chi_X = deg_X * A

    # Step 5: Calculate the middle Betti number b_m = chi(X) - m
    b_middle = chi_X - dim_X

    # Final equation: b_100 = (d1*d2) * (sum) - 100
    print(f"The dimension of the variety is m = {n_ambient} - {len(degrees)} = {dim_X}.")
    print(f"The dimension of the middle cohomology group b_{dim_X} is given by the formula:")
    print(f"b_{dim_X} = chi(X) - {dim_X}")
    print(f"The Euler characteristic chi(X) is computed as:")
    print(f"chi(X) = ({' * '.join(map(str, degrees))}) * [h^{dim_X}] ( (1+h)^{n_ambient+1} / (1+2h)^2 )")
    print(f"The coefficient A = [h^{dim_X}]... is calculated to be {A}.")
    print(f"chi(X) = {deg_X} * {A} = {chi_X}.")
    print(f"Finally, the dimension of the middle cohomology group is:")
    print(f"b_{dim_X} = {chi_X} - {dim_X} = {b_middle}")

solve_cohomology_dimension()