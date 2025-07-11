def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for Y = X1 # X2 # X3.
    """
    print("This script calculates the rank of the kernel K of the Hurewicz homomorphism h*: pi_1(Y) -> H_1(Y).")
    print("The space Y is the connected sum of X1, X2, and X3.")
    print("pi_1(Y) = Z_5 * Z_8 * Z_2, and K is the commutator subgroup of pi_1(Y).\n")
    
    # Orders of the cyclic groups corresponding to X1, X2, X3
    n1 = 5
    n2 = 8
    n3 = 2
    
    print(f"The orders of the fundamental groups pi_1(X1), pi_1(X2), pi_1(X3) are {n1}, {n2}, and {n3} respectively.")
    
    # g_i = |G_i| - 1
    g1 = n1 - 1
    g2 = n2 - 1
    g3 = n3 - 1
    
    print(f"Let g_i = |G_i| - 1. Then we have:")
    print(f"g1 = {n1} - 1 = {g1}")
    print(f"g2 = {n2} - 1 = {g2}")
    print(f"g3 = {n3} - 1 = {g3}\n")
    
    # The rank is calculated using the formula for the commutator subgroup of a free product of n=3 abelian groups:
    # Rank = (2-1)*(g1*g2 + g1*g3 + g2*g3) + (3-1)*(g1*g2*g3)
    
    # Term for pairs (k=2)
    term_k2_summand1 = g1 * g2
    term_k2_summand2 = g1 * g3
    term_k2_summand3 = g2 * g3
    term_k2 = term_k2_summand1 + term_k2_summand2 + term_k2_summand3
    
    print("The formula for the rank involves summing terms for combinations of groups.")
    print("Term for pairs of groups (k=2): (g1*g2 + g1*g3 + g2*g3)")
    print(f"= ({g1}*{g2}) + ({g1}*{g3}) + ({g2}*{g3})")
    print(f"= {term_k2_summand1} + {term_k2_summand2} + {term_k2_summand3} = {term_k2}\n")

    # Term for the triplet (k=3)
    term_k3_product = g1 * g2 * g3
    term_k3 = 2 * term_k3_product
    
    print("Term for the triplet of groups (k=3): 2 * (g1*g2*g3)")
    print(f"= 2 * ({g1}*{g2}*{g3})")
    print(f"= 2 * {term_k3_product} = {term_k3}\n")
    
    # Total rank
    rank = term_k2 + term_k3
    
    print("The total rank is the sum of these terms:")
    print(f"Rank = {term_k2} + {term_k3} = {rank}")
    print("\nSo, the rank of the kernel K as a free group is 95.")

solve_rank_of_kernel()