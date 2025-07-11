def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H based on the reasoning that
    the capacity in the context of the strong product is the independence number.
    """
    # The independence number of G = K_m - C_5 is the clique number of C_5.
    # A 5-cycle (C_5) has a maximum clique size of 2.
    omega_C5 = 2
    alpha_G = omega_C5
    
    print("The graph G is K_m with a C_5 removed.")
    print("Its independence number, alpha(G), is the clique number of C_5.")
    print(f"alpha(G) = {alpha_G}")
    print("-" * 20)

    # The independence number of H = K_n - C_4 is the clique number of C_4.
    # A 4-cycle (C_4) has a maximum clique size of 2.
    omega_C4 = 2
    alpha_H = omega_C4
    
    print("The graph H is K_n with a C_4 removed.")
    print("Its independence number, alpha(H), is the clique number of C_4.")
    print(f"alpha(H) = {alpha_H}")
    print("-" * 20)

    # The Shannon capacity of the strong product G⊠H, interpreted as alpha(G⊠H),
    # is the product of the individual independence numbers.
    shannon_capacity = alpha_G * alpha_H
    
    print("The Shannon capacity of the strong product G⊠H simplifies to alpha(G) * alpha(H).")
    print(f"Shannon Capacity(G⊠H) = {alpha_G} * {alpha_H} = {shannon_capacity}")

solve_shannon_capacity()