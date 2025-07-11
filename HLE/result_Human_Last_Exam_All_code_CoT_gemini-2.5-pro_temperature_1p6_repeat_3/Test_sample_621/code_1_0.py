def find_rotation_index(A, B):
    """
    Finds the rotation index 'i' such that B is a cyclic shift of A.
    B = A[i:] + A[:i]
    
    This algorithm has a time complexity of O(n) and is optimal.
    
    Args:
        A: A list of n unique integers.
        B: A list which is a cyclic shift of A.
        
    Returns:
        The integer index 'i'. Returns -1 if inputs are invalid.
    """
    n = len(A)
    if n != len(B) or n == 0:
        return -1 # Or raise an error for invalid input

    # The first element of B must be A[i].
    first_element_of_B = B[0]
    
    try:
        # Find the index of B[0] in A. This is our rotation index 'i'.
        # The .index() method on a list performs a linear search, which is O(n).
        i = A.index(first_element_of_B)
        return i
    except ValueError:
        # This case should not be reached if B is guaranteed to be a rotation of A.
        return -1

# --- Example Usage ---

# Define the two lists, A and B, where B is a rotated version of A.
A = [4, 5, 6, 7, 0, 1, 2, 3]
B = [0, 1, 2, 3, 4, 5, 6, 7]

# Find the rotation index 'i'.
i = find_rotation_index(A, B)

# Print the results and the final equation.
print(f"Given A: {A}")
print(f"Given B: {B}")

if i != -1:
    print(f"\nThe rotation index i is: {i}")
    
    # Construct the parts of the rotated list for the equation
    part1 = A[i:]
    part2 = A[:i]
    
    # Fulfill the requirement to output each number in the final equation.
    # The default string representation of lists serves this purpose.
    print("\nVerification equation:")
    print(f"{B} = {part1} + {part2}")
else:
    print("Could not find a valid rotation index. B is not a cyclic shift of A.")
