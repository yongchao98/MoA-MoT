import math

def is_prime(n):
    """
    Checks if a number is prime. For this problem, primes are within the 0-11 range.
    The primes up to 11 are 2, 3, 5, 7, 11.
    """
    return n in [2, 3, 5, 7, 11]

def apply_vertical_transform(p):
    """
    Applies the vertical transformation rules to a triplet p = [px, py, pz].
    """
    px, py, pz = p[0], p[1], p[2]
    
    # Check if the z value of the previous triplet is prime
    if is_prime(pz):
        # Apply the 'prime' rule set
        next_x = (px - 3 + py) % 12
        next_y = (py + px) % 12
        next_z = (pz * 2 + px) % 12
    else:
        # Apply the 'not prime' rule set
        next_x = (px + 2 - py) % 12
        next_y = (py * 2 - px) % 12
        next_z = (pz + 3 + pz) % 12
        
    return [next_x, next_y, next_z]

def solve_puzzle():
    """
    Solves the puzzle by calculating the missing elements and their sum.
    """
    # Given triplets from the second row
    m11 = [8, 4, 10]
    m12 = [3, 1, 8]
    
    # The problem description is inconsistent. The most logical approach is to
    # assume the vertical transformation rules are correct and apply them to the
    # known second row to find the unknown third row.
    
    # Calculate the first missing triplet M[2][1] from M[1][1]
    # M[1][1] has z=10, which is not prime.
    m21 = apply_vertical_transform(m11)
    
    # Calculate the second missing triplet M[2][2] from M[1][2]
    # M[1][2] has z=8, which is not prime.
    m22 = apply_vertical_transform(m12)
    
    print(f"The first missing triplet [? ? ?] is determined to be: {m21}")
    print(f"The second missing triplet [? ? ?] is determined to be: {m22}")
    
    # Calculate the sum of the elements of the missing triplets
    missing_elements = m21 + m22
    total_sum = sum(missing_elements)
    
    # Create the equation string
    equation_str = " + ".join(map(str, missing_elements))
    
    print("\nThe sum of the missing elements is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

# Run the solver
solve_puzzle()

# For direct answer extraction
missing_triplet1 = apply_vertical_transform([8, 4, 10])
missing_triplet2 = apply_vertical_transform([3, 1, 8])
final_sum = sum(missing_triplet1) + sum(missing_triplet2)
# print(f"<<<{final_sum}>>>")