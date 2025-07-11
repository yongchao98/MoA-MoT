import math

def combinations(n, k):
    """Computes the binomial coefficient n choose k."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def euler_characteristic_hypersurface(n, d):
    """Computes the Euler characteristic of a smooth hypersurface of degree d in CP^n."""
    chi = 0
    for j in range(n):
        chi += combinations(n + 1, j) * ((-d)**(n - 1 - j))
    chi *= d
    return chi

def solve_pi3_rank():
    """
    Calculates the Betti numbers and explains the steps to find the rank of pi_3(X).
    """
    n = 3  # Ambient space is CP^3
    d = 5  # Hypersurface is quintic

    # Step 1: Compute Betti numbers of X
    chi_X = euler_characteristic_hypersurface(n, d)
    
    # b_k(X) = b_k(CP^{n-1}) for k != n-1.
    # For n=3, X is a surface (dim_C = 2). b_k(X) = b_k(CP^2) for k != 2.
    b0_X = 1
    b1_X = 0
    b3_X = 0 # by Poincare duality b3=b1
    b4_X = 1 # by Poincare duality b4=b0

    # From chi(X) = b0 - b1 + b2 - b3 + b4 = 1 - 0 + b2 - 0 + 1 = 2 + b2
    b2_X = chi_X - 2
    
    print("Step 1: Compute the Betti numbers for the quintic hypersurface X.")
    print(f"The Euler characteristic of X is chi(X) = {chi_X}.")
    print(f"The Betti numbers of X are:")
    print(f"b0(X) = {b0_X}")
    print(f"b1(X) = {b1_X}")
    print(f"b2(X) = {b2_X}")
    print(f"b3(X) = {b3_X}")
    print(f"b4(X) = {b4_X}\n")

    print("Step 2: Relate pi_3(X) to a relative homotopy group.")
    print("From the long exact sequence of homotopy groups for the pair (CP^3, X):")
    print("... -> pi_4(CP^3) -> pi_4(CP^3, X) -> pi_3(X) -> pi_3(CP^3) -> ...")
    print("For CP^n, pi_k(CP^n) = pi_k(S^{2n+1}) for k >= 2.")
    print("For CP^3, pi_3(CP^3) = pi_3(S^7) = 0 and pi_4(CP^3) = pi_4(S^7) = 0.")
    print("The sequence simplifies to: 0 -> pi_4(CP^3, X) -> pi_3(X) -> 0.")
    print("This implies pi_3(X) is isomorphic to pi_4(CP^3, X).\n")

    print("Step 3: Relate the relative homotopy group to a relative homology group.")
    print("The rank of pi_3(X) is equal to the rank of pi_4(CP^3, X).")
    print("By the rational Hurewicz theorem for a 2-connected pair (like (CP^3, X)),")
    print("rank(pi_4(CP^3, X)) = dim_Q H_4(CP^3, X; Q).\n")

    print("Step 4: Compute the dimension of the relative homology group.")
    print("From the long exact sequence of homology for the pair (CP^3, X) with Q coefficients:")
    print("... -> H_4(X; Q) -> H_4(CP^3; Q) -> H_4(CP^3, X; Q) -> H_3(X; Q) -> ...")
    print(f"We know H_3(X; Q) = 0 (since b3 = {b3_X}) and H_4(X; Q) is Q (since b4 = {b4_X}).")
    print(f"Also, H_4(CP^3; Q) is Q.")
    print("The sequence becomes: Q -> Q -> H_4(CP^3, X; Q) -> 0.")
    print("The map H_4(X; Q) -> H_4(CP^3; Q) corresponds to inclusion and is multiplication by the degree d=5.")
    print("This map is an isomorphism on Q. Therefore, its cokernel is 0.")
    print("Thus, H_4(CP^3, X; Q) = 0.\n")

    print("Step 5: Final conclusion.")
    rank_pi3 = 0
    print(f"Combining these steps, we have rank(pi_3(X)) = dim_Q H_4(CP^3, X; Q) = {rank_pi3}.")

solve_pi3_rank()