def is_prime(n):
    """Checks if a number is prime (for numbers up to 11)."""
    return n in [2, 3, 5, 7, 11]

def calculate_triplet(left, top):
    """
    Calculates a new triplet based on the left and top neighbors
    using the Vertical Transformation rules.
    """
    x, y, z = left
    prev_x, prev_y, prev_z = top
    
    if is_prime(prev_z):
        # Rule if previous z is prime
        new_x = (x - 3 + prev_y) % 12
        new_y = (y + prev_x) % 12
        new_z = (z * 2 + prev_x) % 12
    else:
        # Rule if previous z is not prime
        new_x = (x + 2 - prev_y) % 12
        new_y = (y * 2 - prev_x) % 12
        new_z = (z + 3 + prev_z) % 12
        
    return [new_x, new_y, new_z]

def solve_matrix_puzzle():
    """
    Solves the puzzle by calculating the missing elements and their sum.
    """
    # Known triplets from the matrix
    m20 = [7, 2, 9]  # M[2][0]
    m11 = [8, 4, 10] # M[1][1]
    m12 = [3, 1, 8]  # M[1][2]
    
    # Calculate the missing triplet M[2][1]
    m21 = calculate_triplet(m20, m11)
    
    # Calculate the missing triplet M[2][2]
    m22 = calculate_triplet(m21, m12)
    
    # Sum the elements of the missing triplets
    total_sum = sum(m21) + sum(m22)
    
    # Prepare the string for the equation
    equation_str = " + ".join(map(str, m21 + m22))
    
    print(f"The first missing triplet is: {m21}")
    print(f"The second missing triplet is: {m22}")
    print("\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_matrix_puzzle()
<<<39>>>