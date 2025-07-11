def solve_matrix_puzzle():
    """
    Solves the puzzle by calculating the missing triplets and their sum.
    The logic is based on deriving the correct transformation rules from the provided matrix data.
    """

    # Previous triplet in the middle column (M2)
    m2 = [8, 4, 10]

    # --- Step 1: Calculate the middle triplet of row 3 (M3) ---
    # This is a vertical transformation from M2.
    # The z-value of M2 is 10, which is not prime.
    # We use the vertical transformation rule for non-prime z.
    # Rule derivation:
    # nx = (px - py - 2) % 12  -> (8 - 4 - 2) % 12 = 2
    # ny = (py * 2 - px) % 12 -> (4 * 2 - 8) % 12 = 0
    # nz = (pz * 2) % 12      -> (10 * 2) % 12 = 8
    
    px, py, pz = m2[0], m2[1], m2[2]
    
    m3_x = (px - py - 2) % 12
    m3_y = (py * 2 - px) % 12
    m3_z = (pz * 2) % 12
    
    m3 = [m3_x, m3_y, m3_z]

    # --- Step 2: Calculate the right triplet of row 3 (R3) ---
    # This is a horizontal transformation from the newly calculated M3.
    # M3 = [2, 0, 8]. The sum x + y is 2 + 0 = 2, which is <= 10.
    # We use the horizontal transformation rule for x + y <= 10.
    # nx = (x * 2 + y) % 12 -> (2 * 2 + 0) % 12 = 4
    # ny = (y * 3 - 2) % 12 -> (0 * 3 - 2) % 12 = 10
    # nz = (z * 2) % 12      -> (8 * 2) % 12 = 4

    x, y, z = m3[0], m3[1], m3[2]

    r3_x = (x * 2 + y) % 12
    r3_y = (y * 3 - 2) % 12
    r3_z = (z * 2) % 12

    r3 = [r3_x, r3_y, r3_z]
    
    # --- Step 3: Calculate the sum of all missing elements ---
    missing_elements = m3 + r3
    total_sum = sum(missing_elements)

    # Print the final results as required
    print(f"The first missing triplet is: {m3}")
    print(f"The second missing triplet is: {r3}")
    
    equation_parts = [str(n) for n in missing_elements]
    equation_str = " + ".join(equation_parts)
    
    print(f"The final sum is calculated as: {equation_str} = {total_sum}")
    print(f"<<<{total_sum}>>>")

solve_matrix_puzzle()