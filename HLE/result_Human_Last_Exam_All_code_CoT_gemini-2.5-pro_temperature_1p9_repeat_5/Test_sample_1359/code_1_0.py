def solve_puzzle():
    """
    Calculates the missing elements in the matrix based on the provided rules
    and computes their sum.
    """
    # The starting given triplet in the third row.
    m20 = [7, 2, 9]

    # --- Step 1: Calculate the middle triplet of row 3 (M[2][1]) ---
    x, y, z = m20
    
    # Horizontal transformation from M[2][0] to M[2][1]
    if x + y > 10:
        next_x1 = (x * 3 - y) % 12
        next_y1 = (y * 2 + 4) % 12
        next_z1 = (z + x) % 12
    else: # x + y <= 10
        next_x1 = (x * 2 + y) % 12
        next_y1 = (y * 3 - 2) % 12
        next_z1 = (z * 2) % 12
    
    m21 = [next_x1, next_y1, next_z1]

    # --- Step 2: Calculate the right triplet of row 3 (M[2][2]) ---
    x, y, z = m21
    
    # Horizontal transformation from M[2][1] to M[2][2]
    if x + y > 10:
        next_x2 = (x * 3 - y) % 12
        next_y2 = (y * 2 + 4) % 12
        next_z2 = (z + x) % 12
    else: # x + y <= 10
        next_x2 = (x * 2 + y) % 12
        next_y2 = (y * 3 - 2) % 12
        next_z2 = (z * 2) % 12

    m22 = [next_x2, next_y2, next_z2]
    
    # --- Step 3: Sum the missing elements and print the result ---
    
    missing_elements = m21 + m22
    total_sum = sum(missing_elements)

    print(f"The first missing triplet is: {m21}")
    print(f"The second missing triplet is: {m22}")
    
    # Constructing the equation string
    equation_str = " + ".join(map(str, missing_elements))
    
    print(f"\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_puzzle()
<<<24>>>