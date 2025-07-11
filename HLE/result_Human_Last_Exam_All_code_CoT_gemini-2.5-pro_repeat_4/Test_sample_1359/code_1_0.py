def solve_matrix_puzzle():
    """
    This function solves the puzzle by calculating the missing elements and their sum
    based on the interpreted rules.
    """

    # Helper function to check for primality for numbers relevant to the problem
    def is_prime(n):
        # Primes modulo 12 are 2, 3, 5, 7, 11
        return n in {2, 3, 5, 7, 11}

    # --- Given Matrix Data ---
    # Row 1
    m1 = [8, 4, 10]
    r1 = [3, 1, 8]
    # Row 2
    l2 = [7, 2, 9]

    # --- Step 1: Calculate the middle triplet of the third row ---
    # The source is the left triplet of row 3: l2 = [7, 2, 9]
    # The context is the z-value of the middle triplet of row 2: m1[2] = 10
    # Since 10 is not prime, we use the "x + y <= 10" horizontal rule set.
    x_source, y_source, z_source = l2
    
    # Next x = (x * 2 + y) mod 12
    m2_x = (x_source * 2 + y_source) % 12
    # Next y = (y * 3 - 2) mod 12
    m2_y = (y_source * 3 - 2) % 12
    # Next z = (z * 2) mod 12
    m2_z = (z_source * 2) % 12
    
    m2 = [m2_x, m2_y, m2_z]

    # --- Step 2: Calculate the right triplet of the third row ---
    # The source is the newly calculated middle triplet: m2
    # The context is the z-value of the right triplet of row 2: r1[2] = 8
    # Since 8 is not prime, we use the "x + y <= 10" horizontal rule set again.
    x_source, y_source, z_source = m2

    # Calculate x, y, z before applying the cross-dependency for x
    # Next x = (x * 2 + y) mod 12
    r2_x_temp = (x_source * 2 + y_source) % 12
    # Next y = (y * 3 - 2) mod 12
    r2_y = (y_source * 3 - 2) % 12
    # Next z = (z * 2) mod 12
    r2_z = (z_source * 2) % 12

    # Apply the cross-dependency for x: Next x = (calculated x - middle y - 1) mod 12
    # Here, "middle y" is the y-value of the middle triplet of the current row, which is m2_y.
    r2_x = (r2_x_temp - m2_y - 1) % 12

    r2 = [r2_x, r2_y, r2_z]
    
    # --- Step 3: Sum the missing elements and print the result ---
    missing_elements = m2 + r2
    total_sum = sum(missing_elements)
    
    equation_str = " + ".join(map(str, missing_elements))
    
    print(f"The first missing triplet is: {m2}")
    print(f"The second missing triplet is: {r2}")
    print("\nFinal calculation:")
    print(f"{equation_str} = {total_sum}")

    # The final answer format as requested
    print(f"\n<<<{total_sum}>>>")

solve_matrix_puzzle()