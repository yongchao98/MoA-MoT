def is_prime(n):
    """Checks if a number is prime for the context of this puzzle (0-11)."""
    return n in [2, 3, 5, 7, 11]

def calculate_missing_triplets():
    """
    Calculates the missing triplets and their sum based on the given rules.
    """
    # Known triplets needed for the calculation
    m31 = [7, 2, 9]  # Left of the first missing triplet
    m22 = [8, 4, 10] # Above the first missing triplet
    m23 = [3, 1, 8]  # Above the second missing triplet

    # --- Step 1: Calculate the first missing triplet M(3,2) ---
    print("--- Calculating the first missing triplet M(3,2) ---")

    # Horizontal transformation from M(3,1)
    print(f"Applying horizontal transformation to left triplet: {m31}")
    h_x, h_y, h_z = m31
    if h_x + h_y > 10:
        interim_x = (h_x * 3 - h_y) % 12
        interim_y = (h_y * 2 + 4) % 12
        interim_z = (h_z + h_x) % 12
    else: # x + y <= 10
        interim_x = (h_x * 2 + h_y) % 12
        interim_y = (h_y * 3 - 2) % 12
        interim_z = (h_z * 2) % 12
    interim_triplet1 = [interim_x, interim_y, interim_z]
    print(f"Intermediate triplet is: {interim_triplet1}")

    # Vertical transformation using M(2,2)
    print(f"Applying vertical transformation with top triplet: {m22}")
    curr_x, curr_y, curr_z = interim_triplet1
    prev_x, prev_y, prev_z = m22
    if is_prime(prev_z):
        final_x1 = (curr_x - 3 + prev_y) % 12
        final_y1 = (curr_y + prev_x) % 12
        final_z1 = (curr_z * 2 + prev_x) % 12
    else: # prev_z is not prime
        final_x1 = (curr_x + 2 - prev_y) % 12
        final_y1 = (curr_y * 2 - prev_x) % 12
        final_z1 = (curr_z + 3 + prev_z) % 12
    m32 = [final_x1, final_y1, final_z1]
    print(f"The first calculated missing triplet M(3,2) is: {m32}\n")


    # --- Step 2: Calculate the second missing triplet M(3,3) ---
    print("--- Calculating the second missing triplet M(3,3) ---")

    # Horizontal transformation from the newly found M(3,2)
    print(f"Applying horizontal transformation to left triplet: {m32}")
    h_x, h_y, h_z = m32
    if h_x + h_y > 10:
        interim_x = (h_x * 3 - h_y) % 12
        interim_y = (h_y * 2 + 4) % 12
        interim_z = (h_z + h_x) % 12
    else: # x + y <= 10
        interim_x = (h_x * 2 + h_y) % 12
        interim_y = (h_y * 3 - 2) % 12
        interim_z = (h_z * 2) % 12
    interim_triplet2 = [interim_x, interim_y, interim_z]
    print(f"Intermediate triplet is: {interim_triplet2}")

    # Vertical transformation using M(2,3)
    print(f"Applying vertical transformation with top triplet: {m23}")
    curr_x, curr_y, curr_z = interim_triplet2
    prev_x, prev_y, prev_z = m23
    if is_prime(prev_z):
        final_x2 = (curr_x - 3 + prev_y) % 12
        final_y2 = (curr_y + prev_x) % 12
        final_z2 = (curr_z * 2 + prev_x) % 12
    else: # prev_z is not prime
        final_x2 = (curr_x + 2 - prev_y) % 12
        final_y2 = (curr_y * 2 - prev_x) % 12
        final_z2 = (curr_z + 3 + prev_z) % 12
    m33 = [final_x2, final_y2, final_z2]
    print(f"The second calculated missing triplet M(3,3) is: {m33}\n")


    # --- Step 3: Calculate the sum ---
    print("--- Calculating the final sum ---")
    missing_elements = m32 + m33
    total_sum = sum(missing_elements)
    
    # Create the equation string as requested
    equation_str = " + ".join(map(str, missing_elements))
    
    print(f"The missing elements are: {', '.join(map(str, missing_elements))}")
    print(f"The sum is: {equation_str} = {total_sum}")

if __name__ == '__main__':
    calculate_missing_triplets()
    print("<<<20>>>")