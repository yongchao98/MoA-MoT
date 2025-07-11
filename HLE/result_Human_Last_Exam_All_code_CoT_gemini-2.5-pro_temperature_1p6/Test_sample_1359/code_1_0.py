def solve_matrix_puzzle():
    """
    Calculates the missing elements in the matrix based on the provided rules
    and computes their sum.
    """

    # The first triplet in the third row is given.
    t3_1 = [7, 2, 9]

    def apply_horizontal_transform(triplet):
        """
        Applies the core horizontal transformation rules to a given triplet.
        """
        x, y, z = triplet
        
        # Determine which set of rules to use based on x + y
        if x + y > 10:
            next_x = (x * 3 - y) % 12
            next_y = (y * 2 + 4) % 12
            next_z = (z + x) % 12
        else:  # x + y <= 10
            next_x = (x * 2 + y) % 12
            next_y = (y * 3 - 2) % 12
            next_z = (z * 2) % 12
            
        return [next_x, next_y, next_z]

    # Calculate the second triplet of the third row from the first one.
    t3_2 = apply_horizontal_transform(t3_1)

    # Calculate the third triplet of the third row from the second one.
    t3_3 = apply_horizontal_transform(t3_2)

    # The missing elements are all the numbers in the calculated triplets.
    missing_elements = t3_2 + t3_3
    total_sum = sum(missing_elements)

    # Print the results in a clear format.
    print(f"The first missing triplet is: {t3_2}")
    print(f"The second missing triplet is: {t3_3}")
    
    # Construct and print the sum equation as requested.
    equation_str = " + ".join(map(str, missing_elements))
    print(f"\nThe sum of the missing elements is calculated as:")
    print(f"{equation_str} = {total_sum}")

solve_matrix_puzzle()
<<<24>>>