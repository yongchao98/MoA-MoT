import cmath
import math

def solve():
    """
    This function calculates the number of entries in the character table of PSL(2, 7)
    whose absolute value is strictly greater than 1.
    """
    # Define the two complex character values
    c_val1 = complex(-0.5, math.sqrt(7) / 2)
    c_val2 = complex(-0.5, -math.sqrt(7) / 2)

    # The character table of G = PSL(2, 7)
    character_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, c_val1, c_val2],
        [3, -1, 0, 1, c_val2, c_val1],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    # Count entries with absolute value > 1 for each row
    counts_per_row = []
    for row in character_table:
        count = 0
        for entry in row:
            if abs(entry) > 1:
                count += 1
        counts_per_row.append(count)

    # Calculate the total count
    total_count = sum(counts_per_row)

    # Generate the string for the final equation based on non-zero counts per row
    equation_numbers = [str(c) for c in counts_per_row if c > 0]
    equation_str = " + ".join(equation_numbers)
    
    # Print the explanation and the final equation
    print("The count of entries per row with absolute value strictly greater than 1 are:")
    print(f"Row 1: {counts_per_row[0]}")
    print(f"Row 2: {counts_per_row[1]}")
    print(f"Row 3: {counts_per_row[2]}")
    print(f"Row 4: {counts_per_row[3]}")
    print(f"Row 5: {counts_per_row[4]}")
    print(f"Row 6: {counts_per_row[5]}")
    print("\nThe final equation is the sum of the non-zero counts per row:")
    print(f"{equation_str} = {total_count}")
    print(f"\nTotal number of entries with absolute value > 1: {total_count}")

solve()