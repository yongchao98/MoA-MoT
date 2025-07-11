def solve_e8_tori():
    """
    Calculates the number of F_q-rational maximal tori for a group of type E_8.
    """
    
    # 1. Define the parameters for the reductive group G of type E_8.
    group_type = "E_8"
    
    # The rank of the E_8 root system.
    rank_G = 8
    
    # The number of roots in the E_8 root system.
    num_roots_G = 240
    
    # 2. The dimension of a reductive group is the sum of its rank and the number of roots.
    dim_G = rank_G + num_roots_G
    
    # 3. According to a theorem by Steinberg, the number of F_q-rational maximal tori
    # in a connected reductive group G over F_q equals the number of unipotent elements in G(F_q).
    # For a split group like E_8, this number is q^(dim(G) - rank(G)).
    
    # Calculate the exponent of q.
    exponent = dim_G - rank_G
    
    # 4. Print the explanation and the final result step-by-step.
    print("To find the number of F_q-rational maximal tori of a group G of type E_8, we follow these steps:")
    print("\nStep 1: Use Steinberg's Theorem.")
    print("The theorem states that the number of tori is equal to the number of unipotent elements in G(F_q).")
    
    print("\nStep 2: Use the formula for the number of unipotent elements.")
    print("For a split group, this is given by q^(dim(G) - rank(G)).")

    print("\nStep 3: Define the parameters for E_8.")
    print(f"  - The rank of E_8 is {rank_G}.")
    print(f"  - The number of roots for E_8 is {num_roots_G}.")
    print(f"  - The dimension of E_8 is the sum of rank and roots: {rank_G} + {num_roots_G} = {dim_G}.")

    print("\nStep 4: Calculate the final number.")
    print("The number of tori is q raised to the power of (dimension - rank).")
    print(f"The final equation is: q^({dim_G} - {rank_G})")
    print(f"This simplifies to: q^{exponent}")
    
solve_e8_tori()
