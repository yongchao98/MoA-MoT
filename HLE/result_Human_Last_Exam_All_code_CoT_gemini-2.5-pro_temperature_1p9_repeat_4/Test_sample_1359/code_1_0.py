def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True

def solve_matrix_puzzle():
    """
    Calculates the missing elements in the matrix and their sum based on the given rules.
    """
    # Given triplets from the problem statement
    L3 = [7, 2, 9]
    M2 = [8, 4, 10]
    R2 = [3, 1, 8]
    
    # --- Calculate Row 3, Middle Triplet (M3) ---
    # Previous triplet is M2, Left triplet is L3
    prev_m = M2
    left_m = L3
    
    # Check if previous z (M2[z]) is prime
    if is_prime(prev_m[2]): # M2[z] = 10, not prime
        # This block will not be executed
        x_m3 = (left_m[0] - 3 + prev_m[1]) % 12
        y_m3 = (left_m[1] + prev_m[0]) % 12
        z_m3 = (left_m[2] * 2 + prev_m[0]) % 12
    else:
        # Vertical Transformation, z is not prime
        x_m3 = (left_m[0] + 2 - prev_m[1]) % 12
        y_m3 = (left_m[1] * 2 - prev_m[0]) % 12
        z_m3 = (left_m[2] + 3 + prev_m[2]) % 12

    M3 = [x_m3, y_m3, z_m3]

    # --- Calculate Row 3, Right Triplet (R3) ---
    # Previous triplet is R2, Left triplet is the calculated M3
    prev_r = R2
    left_r = M3
    
    # Check if previous z (R2[z]) is prime
    if is_prime(prev_r[2]): # R2[z] = 8, not prime
        # This block will not be executed
        pass
    else:
        # Vertical Transformation, z is not prime
        # Incorporating the cross-dependency for R3[x]: "The middle triplet's y value influences the right triplet's x value"
        # The 'x' in the formula for Next x comes from the middle triplet's y (left_r[1]).
        x_r3 = (left_r[1] + 2 - prev_r[1]) % 12
        y_r3 = (left_r[1] * 2 - prev_r[0]) % 12
        z_r3 = (left_r[2] + 3 + prev_r[2]) % 12
        
    R3 = [x_r3, y_r3, z_r3]

    # Sum all the missing elements
    total_sum = sum(M3) + sum(R3)
    
    # Print the equation as requested
    print("The missing triplets are M3={} and R3={}".format(M3, R3))
    print("The final sum is calculated as follows:")
    print("{} + {} + {} + {} + {} + {} = {}".format(M3[0], M3[1], M3[2], R3[0], R3[1], R3[2], total_sum))
    
solve_matrix_puzzle()
<<<42>>>