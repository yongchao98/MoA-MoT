def solve_puzzle():
    """
    Determines the missing elements in the matrix and calculates their sum
    based on the provided transformation rules.
    """

    # The starting triplet for the third row is given.
    m_2_0 = [7, 2, 9]

    # --- Calculate the first missing triplet M[2][1] from M[2][0] ---
    x, y, z = m_2_0
    
    # Determine which horizontal rule to apply
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else: # x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
        
    m_2_1 = [next_x, next_y, next_z]

    # --- Calculate the second missing triplet M[2][2] from M[2][1] ---
    x, y, z = m_2_1

    # Determine which horizontal rule to apply
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else: # x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12

    m_2_2 = [next_x, next_y, next_z]
    
    # --- Sum the missing elements and print the result ---
    missing_elements = m_2_1 + m_2_2
    total_sum = sum(missing_elements)

    # Print the final equation as requested
    equation_str = " + ".join(map(str, missing_elements))
    print(f"The missing triplets are {m_2_1} and {m_2_2}.")
    print(f"The sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_puzzle()
<<<24>>>