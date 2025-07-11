import math

def solve_matrix_puzzle():
    """
    This script calculates the missing elements in a 3x3 matrix puzzle
    and finds their sum based on a given set of transformation rules.
    """

    # Helper function to check for primality for numbers within the puzzle's scope.
    def is_prime(n):
        if n in [2, 3, 5, 7, 11]:
            return True
        return False

    # Step 1: Define the known triplets from the grid that are needed for calculations.
    L3 = [7, 2, 9]  # Left triplet, Row 3
    M2 = [8, 4, 10] # Middle triplet, Row 2
    R2 = [3, 1, 8]  # Right triplet, Row 2

    # Announce the start of the process
    print("Starting calculation for missing triplets...\n")
    
    # --- Step 2: Calculate the middle triplet of the third row (M3) ---

    # Part A: Horizontal transformation from L3 = [7, 2, 9]
    x_h, y_h, z_h = L3
    if x_h + y_h > 10:
        inter_x = (x_h * 3 - y_h) % 12
        inter_y = (y_h * 2 + 4) % 12
        inter_z = (z_h + x_h) % 12
    else: # x + y <= 10
        inter_x = (x_h * 2 + y_h) % 12
        inter_y = (y_h * 3 - 2) % 12
        inter_z = (z_h * 2) % 12
    
    inter_M3 = [inter_x, inter_y, inter_z]

    # Part B: Vertical transformation on the intermediate result using M2 = [8, 4, 10]
    ix, iy, iz = inter_M3
    prev_x, prev_y, prev_z = M2
    
    if is_prime(prev_z):
        M3x = (ix - 3 + prev_y) % 12
        M3y = (iy + prev_x) % 12
        M3z = (iz * 2 + prev_x) % 12
    else:  # z is not prime
        M3x = (ix + 2 - prev_y) % 12
        M3y = (iy * 2 - prev_x) % 12
        M3z = (iz + 3 + prev_z) % 12
    
    M3 = [M3x, M3y, M3z]
    print(f"Determined middle triplet (M3): {M3}")

    # --- Step 3: Calculate the right triplet of the third row (R3) ---

    # Part A: Horizontal transformation from the newly calculated M3
    x_h, y_h, z_h = M3
    if x_h + y_h > 10:
        inter_x = (x_h * 3 - y_h) % 12
        inter_y = (y_h * 2 + 4) % 12
        inter_z = (z_h + x_h) % 12
    else: # x + y <= 10
        inter_x = (x_h * 2 + y_h) % 12
        inter_y = (y_h * 3 - 2) % 12
        inter_z = (z_h * 2) % 12
        
    inter_R3 = [inter_x, inter_y, inter_z]

    # Part B: Vertical transformation on the intermediate result using R2 = [3, 1, 8]
    ix, iy, iz = inter_R3
    prev_x, prev_y, prev_z = R2
    
    if is_prime(prev_z):
        R3x = (ix - 3 + prev_y) % 12
        R3y = (iy + prev_x) % 12
        R3z = (iz * 2 + prev_x) % 12
    else:  # z is not prime
        R3x = (ix + 2 - prev_y) % 12
        R3y = (iy * 2 - prev_x) % 12
        R3z = (iz + 3 + prev_z) % 12
        
    R3 = [R3x, R3y, R3z]
    print(f"Determined right triplet (R3): {R3}\n")

    # --- Step 4: Sum the missing elements and display the result ---
    
    missing_elements = M3 + R3
    total_sum = sum(missing_elements)
    
    print("The final sum is calculated from the missing elements.")
    # Printing each number in the final equation as requested
    equation_str = " + ".join(map(str, missing_elements))
    print(f"Final Equation: {equation_str} = {total_sum}")

solve_matrix_puzzle()