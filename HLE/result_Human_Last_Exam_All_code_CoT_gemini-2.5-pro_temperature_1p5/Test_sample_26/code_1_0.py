import math

def solve():
    """
    Calculates the rank of the third homotopy group pi_3(X) for a smooth
    quintic hypersurface X in CP^3.
    """
    # Parameters of the problem
    # n is the dimension of the ambient complex projective space CP^n
    n = 3
    # d is the degree of the quintic hypersurface
    d = 5

    # Step 1: Compute Betti numbers of X.
    # b_k(X) are the ranks of the homology groups H_k(X).
    # X is connected.
    b0 = 1
    # By Lefschetz Hyperplane Theorem, X is simply connected, so b1=0.
    b1 = 0
    # For a hypersurface in CP^n, b_k matches CP^{n-1} for k != n-1.
    # Here n=3, dim(X)=2. We compare to CP^2. b_3(X)=b_3(CP^2)=0.
    b3 = 0
    # By Poincare Duality for a real 4-manifold, b4=b0.
    b4 = 1

    # To find b2, we compute the Euler characteristic chi(X).
    # Formula for a smooth hypersurface of degree d in CP^n is:
    # chi = ((1-d)^(n+1) - 1)/d + (n+1)
    chi_X = (pow(1 - d, n + 1) - 1) // d + (n + 1)

    # From the definition of Euler characteristic:
    # chi(X) = b0 - b1 + b2 - b3 + b4
    # So, b2 = chi(X) - (b0 - b1 - b3 + b4)
    b2 = chi_X - (b0 - b1 - b3 + b4)

    # Step 2: Calculate the rank of pi_3(X) using the formula:
    # rank(pi_3(X)) = rank(H_3(X)) + rank(ker(cup: Lambda^2 H^2 -> H^4))
    
    # rank(H_3(X)) is the third Betti number.
    rank_H3 = b3

    # rank(ker(cup)) = dim(Lambda^2 H^2) - dim(im(cup))
    # The cup product to H^4 is surjective, so dim(im(cup)) = b4 = 1.
    dim_im_cup = 1
    
    # dim(Lambda^2 H^2) = C(b2, 2)
    dim_lambda2_H2 = math.comb(b2, 2)

    rank_ker_cup = dim_lambda2_H2 - dim_im_cup
    
    rank_pi3 = rank_H3 + rank_ker_cup

    # Step 3: Print all numbers in the final calculation.
    print(f"The rank of H_3(X) is b3 = {rank_H3}.")
    print(f"The second Betti number of X is b2 = {b2}.")
    print(f"The dimension of Lambda^2 H^2(X, Q) is C({b2}, 2) = {dim_lambda2_H2}.")
    print(f"The dimension of the image of the cup product is {dim_im_cup}.")
    print(f"The rank of the kernel of the cup product is {dim_lambda2_H2} - {dim_im_cup} = {rank_ker_cup}.")
    print("\nThe final equation is rank(pi_3(X)) = rank(H_3(X)) + rank(ker(cup)).")
    print(f"rank(pi_3(X)) = {rank_H3} + {rank_ker_cup} = {rank_pi3}")

solve()