def solve_puzzle():
    """
    Solves the puzzle by calculating the missing triplets and their sum
    based on the provided rules.
    """

    def transform_horizontal(triplet):
        """Applies the horizontal transformation rule to a triplet."""
        x, y, z = triplet
        if x + y > 10:
            next_x = (x * 3 - y) % 12
            next_y = (y * 2 + 4) % 12
            next_z = (z + x) % 12
        else:
            next_x = (x * 2 + y) % 12
            next_y = (y * 3 - 2) % 12
            next_z = (z * 2) % 12
        return [next_x, next_y, next_z]

    # The starting triplet for the third row is given.
    r3c1 = [7, 2, 9]
    print(f"Given left triplet of row 3: {r3c1}")

    # Calculate the middle triplet of row 3 from the left triplet.
    r3c2 = transform_horizontal(r3c1)
    print(f"Calculated middle triplet of row 3: {r3c2}")

    # Calculate the right triplet of row 3 from the middle triplet.
    r3c3 = transform_horizontal(r3c2)
    print(f"Calculated right triplet of row 3: {r3c3}")

    # The missing elements are the numbers in the calculated triplets.
    missing_elements = r3c2 + r3c3
    total_sum = sum(missing_elements)

    # Print the final sum equation.
    equation_str = " + ".join(map(str, missing_elements))
    print(f"Sum of missing elements: {equation_str} = {total_sum}")

    # Final answer in the specified format
    print(f"<<<{total_sum}>>>")

solve_puzzle()