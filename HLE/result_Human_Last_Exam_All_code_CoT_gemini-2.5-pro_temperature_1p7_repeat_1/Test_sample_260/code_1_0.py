def solve():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Orders of the cyclic groups from the fundamental groups of X1, X2, X3
    n1 = 5
    n2 = 8
    n3 = 2
    
    # Order of the first homology group H_1(Y)
    order_H = n1 * n2 * n3
    
    # Euler characteristic of the base CW-complex W
    # 1 0-cell, 3 1-cells, 3 2-cells
    chi_W = 1 - 3 + 3
    
    # Euler characteristic of the covering space hat_W
    chi_hat_W = order_H * chi_W
    
    # Ranks of the kernels of the norm element multiplications
    # rank_ker_Na = order_H * (1 - 1/n1)
    rank_ker_N1 = order_H * (n1 - 1) // n1
    # rank_ker_Nb = order_H * (1 - 1/n2)
    rank_ker_N2 = order_H * (n2 - 1) // n2
    # rank_ker_Nc = order_H * (1 - 1/n3)
    rank_ker_N3 = order_H * (n3 - 1) // n3
    
    # The second Betti number b2(hat_W) is the sum of these ranks
    b2_hat_W = rank_ker_N1 + rank_ker_N2 + rank_ker_N3
    
    # The rank r of the free group K = Ker(h_*) is given by the formula:
    # r = b2(hat_W) - (chi(hat_W) - b0(hat_W))
    # Since hat_W is connected, b0(hat_W) = 1.
    b0_hat_W = 1
    rank = b2_hat_W - (chi_hat_W - b0_hat_W)
    
    print("The final calculation for the rank is:")
    print(f"rank = b2 - (chi_hat - b0) = {b2_hat_W} - ({chi_hat_W} - {b0_hat_W}) = {rank}")

solve()