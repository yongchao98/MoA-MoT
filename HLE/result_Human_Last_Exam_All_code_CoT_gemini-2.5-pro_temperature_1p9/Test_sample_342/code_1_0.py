def solve_sequence():
    """
    This function calculates the first 4 terms of the sequence D_g,
    representing the smallest degree of a finite etale cover of A_g
    needed to make the universal polarization symmetric.
    """

    # For g = 1, 2, 3, the value of D_g is the number of odd theta characteristics,
    # as this corresponds to the size of the smallest orbit of Sp(2g, Z).
    # The number of odd theta characteristics for genus g is given by 2^(g-1) * (2^g - 1).

    # Calculate D_1
    g = 1
    d1 = 2**(g-1) * (2**g - 1)

    # Calculate D_2
    g = 2
    d2 = 2**(g-1) * (2**g - 1)

    # Calculate D_3
    g = 3
    d3 = 2**(g-1) * (2**g - 1)

    # For g = 4, the orbits of Sp(8, Z) on theta characteristics are known.
    # The odd characteristics form a single orbit of size 120.
    # The even characteristics split into two orbits of sizes 1 and 135.
    # D_4 is the minimal orbit size.
    d4 = 1

    print("The sequence of the first 4 terms of D_g is:")
    print(f"D_1 = {d1}")
    print(f"D_2 = {d2}")
    print(f"D_3 = {d3}")
    print(f"D_4 = {d4}")

solve_sequence()