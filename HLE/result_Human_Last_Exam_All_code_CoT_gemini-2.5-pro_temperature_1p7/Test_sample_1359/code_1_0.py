import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

# --- Core Transformation Rules ---

def horizontal_transform(t):
    """Applies horizontal transformation to a triplet."""
    x, y, z = t
    if x + y > 10:
        # Rule: x+y > 10
        nx = (x * 3 - y) % 12
        ny = (y * 2 + 4) % 12
        nz = (z + x) % 12
    else:
        # Rule: x+y <= 10
        nx = (x * 2 + y) % 12
        ny = (y * 3 - 2) % 12
        nz = (z * 2) % 12
    return [nx, ny, nz]

def vertical_transform(t):
    """Applies vertical transformation to a triplet."""
    x, y, z = t
    if is_prime(z):
        # Rule: previous z is prime
        nx = (x - 3 + y) % 12
        ny = (y + x) % 12
        nz = (z * 2 + x) % 12
    else:
        # Rule: previous z is not prime
        nx = (x + 2 - y) % 12
        ny = (y * 2 - x) % 12
        nz = (z + 3 + z) % 12
    return [nx, ny, nz]

def solve():
    """
    Solves the puzzle by calculating the missing elements and their sum.
    """
    # Known Triplets from the problem
    M = [
        [[8, 3, 5], [2, 9, 7], [11, 6, 4]],
        [[5, 7, 2], [8, 4, 10], [3, 1, 8]],
        [[7, 2, 9], None, None]  # Row 2 with missing elements
    ]

    # Calculate M[2][1]
    # Triplet from Above: T_A
    T_A_2_1 = M[1][1] # [8, 4, 10]
    # Triplet from Left: T_L
    T_L_2_1 = M[2][0] # [7, 2, 9]

    # Apply transformations
    T_H_2_1 = horizontal_transform(T_L_2_1) # Result of H on Left
    T_V_2_1 = vertical_transform(T_A_2_1) # Result of V on Above

    # Combine according to deduced model
    nx1 = T_V_2_1[0]
    ny1 = T_H_2_1[1]
    nz1 = (T_H_2_1[2] + T_V_2_1[2]) % 12
    M[2][1] = [nx1, ny1, nz1]

    # Calculate M[2][2]
    # Triplet from Above: T_A
    T_A_2_2 = M[1][2] # [3, 1, 8]
    # Triplet from Left: T_L (this is the one we just calculated)
    T_L_2_2 = M[2][1]

    # Apply transformations
    T_H_2_2 = horizontal_transform(T_L_2_2) # Result of H on Left
    T_V_2_2 = vertical_transform(T_A_2_2) # Result of V on Above

    # Combine according to deduced model
    # Get base x from vertical transform
    nx2_base = T_V_2_2[0]
    # Apply special cross-dependency for column 2
    middle_y = M[2][1][1]
    nx2 = (nx2_base - middle_y) % 12
    
    ny2 = T_H_2_2[1]
    nz2 = (T_H_2_2[2] + T_V_2_2[2]) % 12
    M[2][2] = [nx2, ny2, nz2]
    
    # Extract missing elements
    missing_elements = M[2][1] + M[2][2]
    
    # Calculate sum
    total_sum = sum(missing_elements)

    # Print results
    print(f"The missing triplet at position (2,1) is: {M[2][1]}")
    print(f"The missing triplet at position (2,2) is: {M[2][2]}")
    
    # Print the equation as requested
    equation_str = ' + '.join(map(str, missing_elements))
    print(f"The sum of the missing elements is:")
    print(f"{equation_str} = {total_sum}")


solve()
<<<30>>>