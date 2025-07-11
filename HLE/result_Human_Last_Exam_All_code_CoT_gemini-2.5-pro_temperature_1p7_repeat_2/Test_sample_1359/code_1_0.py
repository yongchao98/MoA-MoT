import math

def is_prime(n):
    """Checks if a number is a prime number (2, 3, 5, 7, 11)."""
    return n in [2, 3, 5, 7, 11]

def horizontal_transform(triplet):
    """Applies the horizontal transformation rule to a triplet."""
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

def vertical_transform(triplet, prev_triplet):
    """Applies the vertical transformation rule to a triplet."""
    x, y, z = triplet
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

def solve_puzzle():
    """
    Solves the puzzle by calculating the missing triplets and their sum.
    """
    # The known part of the matrix
    M = [
        [[8, 3, 5], [2, 9, 7], [11, 6, 4]],
        [[5, 7, 2], [8, 4, 10], [3, 1, 8]],
        [[7, 2, 9], None, None]
    ]

    # Calculate the first missing triplet M[2][1]
    left_triplet_1 = M[2][0]
    top_triplet_1 = M[1][1]
    
    intermediate_triplet_1 = horizontal_transform(left_triplet_1)
    missing_triplet_1 = vertical_transform(intermediate_triplet_1, top_triplet_1)
    
    # Update the matrix with the found triplet
    M[2][1] = missing_triplet_1

    # Calculate the second missing triplet M[2][2]
    left_triplet_2 = M[2][1]
    top_triplet_2 = M[1][2]

    intermediate_triplet_2 = horizontal_transform(left_triplet_2)
    missing_triplet_2 = vertical_transform(intermediate_triplet_2, top_triplet_2)
    
    # Sum the elements of the two missing triplets
    all_missing_elements = missing_triplet_1 + missing_triplet_2
    total_sum = sum(all_missing_elements)
    
    # Print the result in the specified format
    equation_str = " + ".join(map(str, all_missing_elements))
    print(f"{equation_str} = {total_sum}")
    print(f"<<<{total_sum}>>>")

# Run the solver
solve_puzzle()