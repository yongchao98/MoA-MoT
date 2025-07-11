import math

def solve_pi3_rank():
    """
    Calculates the rank of the third homotopy group of a smooth quintic
    hypersurface in CP^3.
    """
    
    # Step 1: Define parameters and compute the Euler characteristic.
    # X is a smooth quintic hypersurface in CP^3.
    # d is the degree of the hypersurface.
    # n is the complex dimension of the hypersurface.
    d = 5
    n = 2

    print("Step 1: Compute the topological invariants of the quintic hypersurface X.")
    # The formula for the Euler characteristic chi of a smooth hypersurface of
    # degree d in CP^3 is given by d * (d^2 - 4*d + 6).
    chi = d * (d**2 - 4*d + 6)
    
    print(f"The degree of the hypersurface is d = {d}.")
    print(f"The Euler characteristic is chi = {d} * ({d}^2 - 4*{d} + 6) = {chi}.")

    # By the Lefschetz hyperplane theorem, X is simply connected, so b1 = 0.
    # By Poincare duality, b3 = b1, so b3 = 0.
    # X is connected, so b0 = 1. The fundamental class gives b4 = 1.
    b0 = 1
    b1 = 0
    b3 = 0
    b4 = 1
    
    # We compute the second Betti number b2 using the Euler-Poincare formula:
    # chi = b0 - b1 + b2 - b3 + b4 = 1 - 0 + b2 - 0 + 1 = b2 + 2
    b2 = chi - 2
    print(f"The Betti numbers of X are: b0={b0}, b1={b1}, b3={b3}, b4={b4}.")
    print(f"From the Euler characteristic, we find the second Betti number: b2 = chi - 2 = {chi} - 2 = {b2}.")
    print("-" * 20)

    # Step 2-4: Use rational homotopy theory to find the rank of pi_3(X).
    print("Step 2: Apply results from rational homotopy theory.")
    # X is a simply connected compact Kähler manifold, hence it is a "formal space".
    # For a formal space with H3(X)=0, the rank of pi_3(X) is given by the
    # dimension of the kernel of the cup product map:
    # cup: Lambda^2 H^2(X, Q) -> H^4(X, Q)
    # The dimensions of these spaces are determined by the Betti numbers.
    dim_H2 = b2
    dim_H4 = b4

    # The dimension of the domain Lambda^2 H^2(X, Q) is C(b2, 2).
    dim_Lambda2_H2 = math.comb(dim_H2, 2)
    print(f"The dimension of H^2(X, Q) is b2 = {dim_H2}.")
    print(f"The dimension of the domain of the cup product map, Lambda^2 H^2(X, Q), is C({dim_H2}, 2) = {dim_Lambda2_H2}.")

    # For a Kähler manifold, the intersection form (Poincare dual to cup product)
    # is non-degenerate. This implies the cup product map is surjective onto H^4(X, Q).
    dim_image = dim_H4
    print(f"The cup product map is surjective, so its image has dimension b4 = {dim_image}.")

    # By the rank-nullity theorem, dim(ker) = dim(domain) - dim(image).
    rank_pi3 = dim_Lambda2_H2 - dim_image

    print("\nStep 3: Calculate the rank of pi_3(X).")
    print("The rank is the dimension of the kernel of the cup product map.")
    # Final equation format requested by the user
    print(f"rank(pi_3(X)) = dim(Lambda^2 H^2) - dim(H^4) = C({dim_H2}, 2) - {dim_image}")
    print(f"The final equation is: {dim_Lambda2_H2} - {dim_image} = {rank_pi3}")


solve_pi3_rank()