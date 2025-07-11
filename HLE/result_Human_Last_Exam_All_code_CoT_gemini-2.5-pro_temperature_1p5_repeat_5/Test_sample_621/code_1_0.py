def find_rotation_index(A, B):
    """
    Given two lists A and B, where B is a rotation of A, this function
    finds the rotation index i such that B = A[i:] + A[:i].

    The algorithm has a time complexity of O(n) and is asymptotically optimal.
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same length.")
        return

    # The core of the algorithm: find the index of B's first element in A.
    # This takes O(n) time.
    try:
        first_element_of_B = B[0]
        i = A.index(first_element_of_B)
    except ValueError:
        print(f"Error: The first element of B ({B[0]}) was not found in A.")
        print("This violates the problem assumption that B is a rotation of A.")
        return

    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"The rotation index is i = {i}")
    
    # Reconstruct the equation B = A[i:] + A[:i] to demonstrate the result.
    part1 = A[i:]
    part2 = A[:i]
    
    print("\nDemonstrating the equation B = A[i:] + A[:i]:")
    print(f"A[i:] where i={i} is {part1}")
    print(f"A[:i] where i={i} is {part2}")
    
    # Final equation outputting each number
    print(f"\nFinal Equation:")
    print(f"{B} = {part1} + {part2}")


# Example usage:
# A = [10, 20, 30, 40, 50, 60]
# B = [40, 50, 60, 10, 20, 30] -> i=3
A_example = [10, 20, 30, 40, 50, 60]
B_example = [40, 50, 60, 10, 20, 30]
find_rotation_index(A_example, B_example)
