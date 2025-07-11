def is_prime(n):
    """Checks if a number is prime."""
    if n < 2:
        return False
    primes = {2, 3, 5, 7, 11}
    return n in primes

def calculate_missing_triplets():
    """
    Calculates the missing triplets and their sum based on the given rules.
    """
    # Given matrix values
    L3 = [7, 2, 9]
    M2 = [8, 4, 10]
    R2 = [3, 1, 8]

    # --- Step 1: Calculate the middle triplet of Row 3 (M3) ---
    print("Calculating the middle triplet of Row 3 (M3)...")
    
    # Horizontal transformation from L3
    x, y, z = L3
    print(f"Applying horizontal transformation on L3 = {L3}")
    if x + y > 10:
        temp_x = (x * 3 - y) % 12
        temp_y = (y * 2 + 4) % 12
        temp_z = (z + x) % 12
    else: # x + y <= 10
        temp_x = (x * 2 + y) % 12
        temp_y = (y * 3 - 2) % 12
        temp_z = (z * 2) % 12
    temp_M3 = [temp_x, temp_y, temp_z]
    print(f"Resulting temporary triplet (temp_M3) = {temp_M3}")

    # Vertical transformation using M2
    x, y, z = temp_M3
    prev_x, prev_y, prev_z = M2
    print(f"Applying vertical transformation on temp_M3 using M2 = {M2}")
    if is_prime(prev_z):
        final_x_m3 = (x - 3 + prev_y) % 12
        final_y_m3 = (y + prev_x) % 12
        final_z_m3 = (z * 2 + prev_x) % 12
    else: # prev_z is not prime
        final_x_m3 = (x + 2 - prev_y) % 12
        final_y_m3 = (y * 2 - prev_x) % 12
        final_z_m3 = (z + 3 + prev_z) % 12
    M3 = [final_x_m3, final_y_m3, final_z_m3]
    print(f"Final middle triplet M3 = {M3}\n")

    # --- Step 2: Calculate the right triplet of Row 3 (R3) ---
    print("Calculating the right triplet of Row 3 (R3)...")
    
    # Horizontal transformation from M3
    x, y, z = M3
    print(f"Applying horizontal transformation on M3 = {M3}")
    if x + y > 10:
        temp_x = (x * 3 - y) % 12
        temp_y = (y * 2 + 4) % 12
        temp_z = (z + x) % 12
    else: # x + y <= 10
        temp_x = (x * 2 + y) % 12
        temp_y = (y * 3 - 2) % 12
        temp_z = (z * 2) % 12
    temp_R3 = [temp_x, temp_y, temp_z]
    print(f"Resulting temporary triplet (temp_R3) = {temp_R3}")

    # Vertical transformation using R2
    x, y, z = temp_R3
    prev_x, prev_y, prev_z = R2
    print(f"Applying vertical transformation on temp_R3 using R2 = {R2}")
    if is_prime(prev_z):
        final_x_r3 = (x - 3 + prev_y) % 12
        final_y_r3 = (y + prev_x) % 12
        final_z_r3 = (z * 2 + prev_x) % 12
    else: # prev_z is not prime
        final_x_r3 = (x + 2 - prev_y) % 12
        final_y_r3 = (y * 2 - prev_x) % 12
        final_z_r3 = (z + 3 + prev_z) % 12
    R3 = [final_x_r3, final_y_r3, final_z_r3]
    print(f"Final right triplet R3 = {R3}\n")

    # --- Step 3: Calculate and print the sum ---
    total_sum = sum(M3) + sum(R3)
    
    print("The missing elements are:")
    print(f"M3: {M3[0]}, {M3[1]}, {M3[2]}")
    print(f"R3: {R3[0]}, {R3[1]}, {R3[2]}")
    
    print("\nCalculating the sum of the missing elements:")
    equation_parts = [str(n) for n in M3] + [str(n) for n in R3]
    print(f"{' + '.join(equation_parts)} = {total_sum}")

if __name__ == "__main__":
    calculate_missing_triplets()
    print("\n<<<20>>>")
