import math

def solve_homotopy_rank():
    """
    Calculates the rank of the third homotopy group of a smooth quintic hypersurface in CP^3.
    """
    
    # 1. Properties of the hypersurface X
    d = 5  # degree of the hypersurface
    n_ambient = 3  # dimension of the ambient projective space CP^3
    dim_X = n_ambient - 1 # dimension of X
    
    print("Let X be a smooth quintic hypersurface in CP^3.")
    print(f"The complex dimension of X is {dim_X}, so its real dimension is {2 * dim_X}.")
    
    # 2. Basic topological invariants from Lefschetz Hyperplane Theorem and Poincaré Duality
    print("\nStep 1: Determine basic topological invariants of X.")
    print("By the Lefschetz Hyperplane Theorem, pi_1(X) is isomorphic to pi_1(CP^3), which is 0.")
    print("Since X is simply connected (pi_1(X) = 0), its first Betti number b_1 is 0.")
    print("By Poincaré Duality for a 4-manifold, b_3 = b_(4-3) = b_1, so b_3 = 0.")
    print("This implies that the homology group H_3(X, Z) = 0.")

    b1 = 0
    b3 = 0
    
    # 3. Calculate the second Betti number, b_2
    print("\nStep 2: Calculate the second Betti number b_2(X).")
    
    # K_X^2 = d
    K_X_sq = d
    print(f"The square of the canonical class of X is K_X^2 = d = {K_X_sq}.")

    # p_g = h^{2,0} = C(d-1, n)
    p_g = math.comb(d - 1, n_ambient)
    # q = h^{1,0} = 0
    q = 0
    # chi_O_X = p_g - q + 1
    chi_O_X = p_g - q + 1
    print(f"The geometric genus is p_g = C({d-1}, {n_ambient}) = {p_g}.")
    print(f"The arithmetic genus is chi(O_X) = p_g - q + 1 = {p_g} - {q} + 1 = {chi_O_X}.")

    # Noether's formula: 12 * chi_O_X = K_X^2 + chi_top
    chi_top = 12 * chi_O_X - K_X_sq
    print(f"By Noether's formula, 12 * chi(O_X) = K_X^2 + chi_top(X).")
    print(f"So, chi_top(X) = 12 * {chi_O_X} - {K_X_sq} = {chi_top}.")

    # chi_top = b0 - b1 + b2 - b3 + b4
    # b0=1, b4=1
    b0 = 1
    b4 = 1
    b2 = chi_top - (b0 - b1 - b3 + b4)
    print(f"The Euler characteristic is also given by the sum of Betti numbers:")
    print(f"chi_top(X) = b0 - b1 + b2 - b3 + b4")
    print(f"{chi_top} = {b0} - {b1} + b2 - {b3} + {b4} = {b0 - b1 - b3 + b4} + b2")
    print(f"This gives b2 = {chi_top} - {b0 - b1 - b3 + b4} = {b2}.")
    
    # 4. Connect to homotopy groups via Hurewicz Theorem
    print("\nStep 3: Relate homology to homotopy.")
    k = b2
    print(f"Since X is 1-connected, the Hurewicz Theorem states that pi_2(X) is isomorphic to H_2(X).")
    print(f"H_2(X) is a free abelian group of rank b2 = {k}.")
    print(f"So, pi_2(X) is isomorphic to Z^{k}.")
    
    # 5. Use Whitehead's theorem for pi_3
    print("\nStep 4: Compute the rank of pi_3(X).")
    print("For a 1-connected 4-manifold, a theorem by J.H.C. Whitehead states:")
    print("pi_3(X) is isomorphic to H_3(X) + Gamma(pi_2(X)).")
    print(f"We have H_3(X) = 0 and pi_2(X) = Z^{k}, with k = {k}.")
    print(f"The rank of the free part of Gamma(Z^k) is given by the binomial coefficient C(k, 2).")
    
    # 6. Final Calculation
    rank = math.comb(k, 2)
    print("\nFinal Calculation:")
    print(f"Rank(pi_3(X)) = C({k}, 2) = ({k} * ({k}-1)) / 2 = {rank}")

solve_homotopy_rank()