import cmath

def solve():
    """
    Calculates the number of entries in the character table of PSL(2,7)
    with an absolute value strictly greater than 1.
    """

    # Define the complex numbers alpha and beta found in the character table
    alpha = (-1 + 1j * cmath.sqrt(7)) / 2
    beta = (-1 - 1j * cmath.sqrt(7)) / 2

    # The character table of G = PSL(2,7)
    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, alpha, beta],
        [3, -1, 0, 1, beta, alpha],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    counts_per_character = []
    total_count = 0

    # Iterate through each character (row) in the table
    for char_row in char_table:
        row_count = 0
        # Iterate through each entry in the row
        for entry in char_row:
            if abs(entry) > 1:
                row_count += 1
        counts_per_character.append(row_count)
    
    total_count = sum(counts_per_character)
    
    # Format the final equation showing the sum of counts per character
    equation_str = " + ".join(map(str, counts_per_character))
    
    print(f"The number of entries in each row with absolute value greater than 1 are: {counts_per_character}")
    print("The final count is the sum of these numbers:")
    print(f"{equation_str} = {total_count}")

solve()