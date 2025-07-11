import math

def solve_pi3_rank():
    """
    Calculates the rank of the third homotopy group of a smooth quintic hypersurface in CP^3.
    """
    # Step 1: Define parameters for the smooth quintic hypersurface in CP^3
    n = 3  # Dimension of the ambient projective space CP^n
    d = 5  # Degree of the hypersurface

    print(f"The hypersurface X is a smooth surface of degree d={d} in CP^{n}.")
    print(f"The complex dimension of X is k = n-1 = {n-1}.")
    print(f"The real dimension of X is 2k = {2*(n-1)}.\n")

    # Step 2: Calculate the Euler characteristic chi(X)
    # For a hypersurface of degree d in CP^n, the Euler characteristic can be computed
    # from its Chern classes. For a quintic in CP^3, this is a known result.
    # The coefficient of h^2 in the expansion of (1+h)^4 / (1+5h) is 11.
    # chi(X) = coefficient * d
    c2_coeff = 11
    chi = c2_coeff * d
    print(f"The Euler characteristic chi(X) is calculated from Chern classes.")
    print(f"For a quintic in CP^3, chi(X) = {c2_coeff} * d = {c2_coeff} * {d} = {chi}.\n")

    # Step 3: Calculate the Betti numbers of X
    # b_i = rank(H_i(X, Z))
    b0 = 1
    b1 = 0  # By Lefschetz hyperplane theorem
    b4 = 1  # By Poincare duality, b_4 = b_0
    b3 = 0  # By Poincare duality, b_3 = b_1
    
    # The Euler characteristic is the alternating sum of Betti numbers.
    # chi = b0 - b1 + b2 - b3 + b4
    # chi = 1 - 0 + b2 - 0 + 1 = b2 + 2
    b2 = chi - (b0 - b1 + b3 + b4)
    print("The Betti numbers of X are:")
    print(f"b0 = {b0}")
    print(f"b1 = {b1}")
    print(f"b3 = {b3}")
    print(f"b4 = {b4}")
    print(f"The second Betti number b2 is derived from the Euler characteristic:")
    print(f"b2 = chi(X) - (b0 - b1 + b3 + b4) = {chi} - ({b0} - {b1} + {b3} + {b4}) = {b2}\n")

    # Step 4: Apply the theorem for the rank of pi_3(X)
    # For a simply-connected algebraic surface X of general type, the rank of pi_3(X)
    # is the rank of the Whitehead group Gamma(H_2(X)).
    # Since H_2(X) is isomorphic to Z^b2, this rank is the binomial coefficient C(b2, 2).
    k = b2
    if k < 2:
        rank_pi3 = 0
    else:
        rank_pi3 = math.comb(k, 2)
    
    print("For a simply-connected algebraic surface like X (which is of general type),")
    print("the rank of the third homotopy group pi_3(X) is given by the binomial coefficient 'b2 choose 2'.")
    
    numerator = k * (k - 1)
    
    print(f"Rank(pi_3(X)) = C({k}, 2) = ({k} * ({k}-1)) / 2 = {numerator} / 2 = {rank_pi3}")

solve_pi3_rank()