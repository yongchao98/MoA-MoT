from math import comb

def solve_pi3_rank():
    """
    Calculates the rank of the third homotopy group of a smooth quintic hypersurface in CP^3.
    """
    d = 5 # Degree of the hypersurface
    m = 3 # Dimension of the ambient complex projective space

    print(f"The task is to find the rank of pi_3(X) for a smooth quintic (d={d}) hypersurface X in CP^{m}.")
    print("According to a theorem by Libgober and Wood, the rank is given by the formula:")
    print(f"Rank(pi_{m}(X)) = b_{m-1}(X) - b_{m-1}(CP^{m-1})")
    print(f"For m={m}, this is: Rank(pi_3(X)) = b_2(X) - b_2(CP^2)\n")

    # Step 1: Calculate b_2(X)
    print("First, we calculate the second Betti number, b_2(X).")
    
    # Calculate geometric genus h^{2,0}(X)
    h20 = comb(d - 1, m)
    print(f"The geometric genus h^2,0(X) = C(d-1, m) = C({d-1}, {m}) = {h20}")

    # Calculate chi(O_X)
    chi_OX = 1 + h20
    print(f"The Euler characteristic of the structure sheaf chi(O_X) = 1 + h^2,0(X) = 1 + {h20} = {chi_OX}")

    # Calculate K_X^2
    K_X_sq = (d - m - 1)**2 * d
    print(f"The self-intersection of the canonical class K_X^2 = (d - m - 1)^2 * d = ({d} - {m} - 1)^2 * {d} = {K_X_sq}")
    
    # Calculate topological Euler characteristic chi(X) using Noether's formula
    chi_X = 12 * chi_OX - K_X_sq
    print(f"Using Noether's formula, the topological Euler characteristic chi(X) = 12 * chi(O_X) - K_X^2 = 12 * {chi_OX} - {K_X_sq} = {chi_X}")
    
    # Calculate b_2(X)
    # For a simply connected surface, chi(X) = b_0 - b_1 + b_2 - b_3 + b_4 = 1 - 0 + b_2 - 0 + 1 = b_2 + 2
    b2_X = chi_X - 2
    print(f"For a smooth hypersurface, b_0=1, b_1=0, b_3=0, b_4=1. From chi(X) = b_2(X) + 2, we get b_2(X) = chi(X) - 2 = {chi_X} - 2 = {b2_X}\n")

    # Step 2: Get b_2(CP^2)
    b2_CP2 = 1
    print(f"The second Betti number of CP^2 is b_2(CP^2) = {b2_CP2}.\n")

    # Step 3: Calculate the final rank
    rank = b2_X - b2_CP2
    print("Finally, we compute the rank of pi_3(X):")
    print(f"Rank(pi_3(X)) = b_2(X) - b_2(CP^2) = {b2_X} - {b2_CP2} = {rank}")

solve_pi3_rank()