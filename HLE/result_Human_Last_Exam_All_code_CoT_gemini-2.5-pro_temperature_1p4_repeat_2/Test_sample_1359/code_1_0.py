def is_prime(n):
    """Checks if a number is prime (for values 0-11)."""
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(triplet):
    """Applies the horizontal transformation rule."""
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

def vertical_transform(base_triplet, prev_triplet):
    """Applies the vertical transformation rule."""
    x, y, z = base_triplet
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

def combine(base_triplet, v_triplet, prev_triplet):
    """Combines base and v-transformed triplets based on the meta-rule."""
    b_x, b_y, b_z = base_triplet
    v_x, v_y, v_z = v_triplet
    
    if is_prime(prev_triplet[2]):
        final_x = b_x % 12
        final_y = (b_y - 2) % 12
        final_z = (b_z + 3) % 12
    else:
        final_x = (v_x - 1) % 12
        final_y = v_y % 12
        final_z = (b_z + 2) % 12
    return [final_x, final_y, final_z]

def solve_puzzle():
    """Main function to solve the puzzle."""
    # Given matrix values
    M22 = [8, 4, 10]
    M23 = [3, 1, 8]
    M31 = [7, 2, 9]

    # --- Calculate M32 ---
    # Left: M31, Top: M22
    base_32 = horizontal_transform(M31)
    v_32 = vertical_transform(base_32, M22)
    M32 = combine(base_32, v_32, M22)
    
    # --- Calculate M33 ---
    # Left: M32, Top: M23
    base_33 = horizontal_transform(M32)
    v_33 = vertical_transform(base_33, M23)
    M33 = combine(base_33, v_33, M23)

    # --- Calculate and print the sum ---
    missing_elements = M32 + M33
    total_sum = sum(missing_elements)
    
    equation = " + ".join(map(str, missing_elements))
    print(f"The missing triplets are M32={M32} and M33={M33}.")
    print(f"The sum of the missing elements is calculated as follows:")
    print(f"{equation} = {total_sum}")

solve_puzzle()
<<<22>>>