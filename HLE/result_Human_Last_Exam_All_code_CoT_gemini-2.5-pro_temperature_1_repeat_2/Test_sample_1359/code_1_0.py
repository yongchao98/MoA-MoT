def solve_matrix_puzzle():
    """
    This function calculates the missing elements in the matrix based on the provided rules
    and then prints their sum.
    """

    # The first triplet of the third row is given.
    m31 = [7, 2, 9]

    # --- Calculate the second triplet of the third row: M(3, 2) ---
    x1, y1, z1 = m31
    m32 = [0, 0, 0]

    # Check condition for horizontal transformation
    if (x1 + y1) > 10:
        m32[0] = (x1 * 3 - y1) % 12
        m32[1] = (y1 * 2 + 4) % 12
        m32[2] = (z1 + x1) % 12
    else:
        m32[0] = (x1 * 2 + y1) % 12
        m32[1] = (y1 * 3 - 2) % 12
        m32[2] = (z1 * 2) % 12

    # --- Calculate the third triplet of the third row: M(3, 3) ---
    x2, y2, z2 = m32
    m33 = [0, 0, 0]
    
    # Check condition for horizontal transformation
    if (x2 + y2) > 10:
        x3_base = (x2 * 3 - y2) % 12
        m33[1] = (y2 * 2 + 4) % 12
        m33[2] = (z2 + x2) % 12
    else:
        x3_base = (x2 * 2 + y2) % 12
        m33[1] = (y2 * 3 - 2) % 12
        m33[2] = (z2 * 2) % 12

    # Apply the cross-dependency rule for x3
    # From row 2 analysis: when middle y is 4, the adjustment for x is -5.
    # M(2,2)=[8,4,10], x_base=(8*3-4)%12=8, x_actual=3. Adj = 3-8 = -5.
    # Here, y2 (middle y) is 4, so we use the same adjustment.
    adjustment = -5
    m33[0] = (x3_base + adjustment) % 12
    
    # --- Calculate the sum of the missing elements ---
    missing_elements = m32 + m33
    total_sum = sum(missing_elements)

    # Print the results as requested
    print(f"The second missing triplet is: {m32}")
    print(f"The third missing triplet is: {m33}")
    
    equation_str = " + ".join(map(str, missing_elements))
    print(f"\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_matrix_puzzle()
<<<31>>>