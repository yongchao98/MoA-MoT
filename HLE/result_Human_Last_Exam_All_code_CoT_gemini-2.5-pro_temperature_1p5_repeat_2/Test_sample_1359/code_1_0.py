def is_prime(n):
    """Checks if a number in the context of modulo 12 is prime."""
    # Primes up to 11 are 2, 3, 5, 7, 11
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(triplet):
    """
    Applies the Horizontal Transformation rules.
    Returns a new triplet [x, y, z].
    """
    x, y, z = triplet
    if x + y > 10:
        next_x = (x * 3 - y) % 12
        next_y = (y * 2 + 4) % 12
        next_z = (z + x) % 12
    else: # x + y <= 10
        next_x = (x * 2 + y) % 12
        next_y = (y * 3 - 2) % 12
        next_z = (z * 2) % 12
    return [next_x, next_y, next_z]

def vertical_transform_y(previous_triplet):
    """
    Applies the Vertical Transformation rule for the 'y' component.
    Uses the "rewritten" interpretation where variables on the RHS are from the 'previous' triplet.
    """
    prev_x, prev_y, prev_z = previous_triplet
    if is_prime(prev_z):
        # Next y = (y + previous x) mod 12  ->  (prev_y + prev_x) mod 12
        next_y = (prev_y + prev_x) % 12
    else: # previous z is not prime
        # Next y = (y * 2 - previous x) mod 12 -> (prev_y * 2 - prev_x) mod 12
        next_y = (prev_y * 2 - prev_x) % 12
    # Python's % can give negative results, so we ensure it's in the 0-11 range
    return (next_y + 12) % 12

def solve():
    """
    Calculates the missing elements and their sum based on the hybrid model.
    """
    # Given matrix values needed for calculation
    M31 = [7, 2, 9]
    M22 = [8, 4, 10]
    M23 = [3, 1, 8]

    # --- Calculate M[3][2] ---
    
    # x and z components from M[3][1] horizontally
    h_transform_M31 = horizontal_transform(M31)
    x32 = h_transform_M31[0]
    z32 = h_transform_M31[2]
    
    # y component from M[2][2] vertically
    y32 = vertical_transform_y(M22)
    
    M32 = [x32, y32, z32]

    # --- Calculate M[3][3] ---

    # x and z components from the newly calculated M[3][2] horizontally
    h_transform_M32 = horizontal_transform(M32)
    x33 = h_transform_M32[0]
    z33 = h_transform_M32[2]

    # y component from M[2][3] vertically
    y33 = vertical_transform_y(M23)
    
    M33 = [x33, y33, z33]
    
    missing_elements = M32 + M33
    total_sum = sum(missing_elements)

    print("The missing triplet [? ? ?] in the middle of the last row is:", M32)
    print("The missing triplet [? ? ?] at the right of the last row is:", M33)
    
    # Create the equation string
    equation_str = " + ".join(map(str, missing_elements))
    
    print(f"The sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve()
<<<29>>>