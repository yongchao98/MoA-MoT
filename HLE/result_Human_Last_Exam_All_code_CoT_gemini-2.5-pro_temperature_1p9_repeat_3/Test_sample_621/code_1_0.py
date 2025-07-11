import random

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a rotation of A.
    This is based on the logic that B[0] must equal A[i].

    Args:
        A: The original list of n unique integers.
        B: A list that is a rotation of A.

    Returns:
        The rotation index i, or -1 if B is not a rotation of A
        or lists are empty.
    """
    if not A or not B or len(A) != len(B):
        return -1

    # The first element of B must be the i-th element of A.
    first_element_of_B = B[0]

    # Find the index of this element in A. This is an O(n) operation
    # for a list, which is the bottleneck and makes the algorithm O(n).
    try:
        # Since A contains unique integers, this index is unique.
        i = A.index(first_element_of_B)
        return i
    except ValueError:
        # This case should not be reached if B is guaranteed to be a rotation of A.
        return -1

def main():
    """
    Main function to demonstrate finding the rotation index.
    """
    # Let's create an example.
    n = 10
    A = random.sample(range(1, 101), n) # A list of 10 unique integers.
    
    # Let's pick a random rotation index i.
    i = random.randint(0, n - 1)
    
    # Create the rotated list B.
    B = A[i:] + A[:i]
    
    print(f"Original list A: {A}")
    print(f"List B (rotated from A by i={i}): {B}")
    print("-" * 20)
    
    # Find the rotation index using our algorithm.
    found_i = find_rotation_index(A, B)
    
    if found_i != -1:
        print(f"The algorithm found the rotation index i = {found_i}")
        
        # Verify and display the result as an equation with all numbers.
        print("\nVerifying the result:")
        print(f"Equation: B = A[i:] + A[:i]")
        print(f"B: {B}")
        
        A_part1 = A[found_i:]
        A_part2 = A[:found_i]
        
        print(f"A[{found_i}:] + A[:{found_i}] is: {A_part1} + {A_part2}")
        print(f"Result of concatenation: {A_part1 + A_part2}")
        
        if (A_part1 + A_part2) == B:
            print("\nVerification successful: The lists match.")
        else:
            print("\nVerification failed: The lists do not match.")
            
    else:
        print("Could not find the rotation index.")

if __name__ == "__main__":
    main()
