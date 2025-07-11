def is_prime(n):
    """Checks if a number is prime for the context of this problem (0-11)."""
    # Primes up to 11 are 2, 3, 5, 7, 11
    return n in [2, 3, 5, 7, 11]

def apply_horizontal(triplet):
    """Applies horizontal transformation rules to a triplet."""
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

def apply_vertical(input_triplet, prev_triplet):
    """
    Applies vertical transformation rules.
    - input_triplet is the [x, y, z] to be transformed.
    - prev_triplet is the [x, y, z] from the row above.
    """
    x, y, z = input_triplet
    prev_x, prev_y, prev_z = prev_triplet

    if is_prime(prev_z):
        next_x = (x - 3 + prev_y) % 12
        next_y = (y + prev_x) % 12
        next_z = (z * 2 + prev_x) % 12
    else:
        next_x = (x + 2 - prev_y) % 12
        next_y = (y * 2 - prev_x) % 12
        next_z = (z + 3 + prev_z) % 12
    return [next_x, next_y, next_z]

def solve_matrix_puzzle():
    """
    Solves the puzzle by calculating the missing triplets following a systematic plan.
    """
    # Given known triplets from the matrix needed for the calculation
    T_mid_mid = [8, 4, 10]
    T_mid_right = [3, 1, 8]
    T_bot_left = [7, 2, 9]

    print("Step 1: Calculate the Bottom-Middle triplet at position (2, 1)")
    
    # Horizontal step
    H_intermediate = apply_horizontal(T_bot_left)
    
    # Vertical step
    H_final = apply_vertical(H_intermediate, T_mid_mid)
    
    print(f"The first missing triplet is: {H_final}")
    
    print("\nStep 2: Calculate the Bottom-Right triplet at position (2, 2)")
    
    # Horizontal step
    I_intermediate = apply_horizontal(H_final)

    # Vertical step
    I_final = apply_vertical(I_intermediate, T_mid_right)

    print(f"The second missing triplet is: {I_final}")

    # Combine the elements for the final sum
    missing_elements = H_final + I_final
    total_sum = sum(missing_elements)

    print("\n--- Calculation Summary ---")
    
    # Format the sum equation as requested
    equation_str = " + ".join(map(str, missing_elements))
    print(f"The sum of the missing elements is: {equation_str} = {total_sum}")


# Run the solver
solve_matrix_puzzle()
<<<20>>>