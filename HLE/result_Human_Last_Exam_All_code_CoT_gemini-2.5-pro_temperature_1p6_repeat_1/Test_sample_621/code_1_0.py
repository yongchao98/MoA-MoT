def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    The relationship is B = A[i:] + A[:i].

    Args:
        A: A list with n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i, or None if inputs are invalid.
    """
    # Basic validation for non-empty lists of the same size.
    if not A or not B or len(A) != len(B):
        print("Error: Input lists must be non-empty and have the same length.")
        return None

    n = len(A)
    if n == 0:
        print("Rotation index is 0 for empty lists.")
        return 0
        
    # The first element of B must correspond to A[i].
    target_element = B[0]

    # Find the index 'i' of the target_element in list A.
    # This list.index() method performs a linear scan, which is an O(n) operation.
    try:
        i = A.index(target_element)
    except ValueError:
        # This case should not happen based on the problem description,
        # but it's good practice for robust code.
        print(f"Error: List B does not seem to be a rotation of A, because {target_element} is not in A.")
        return None

    # Print the results as requested.
    print(f"Given A = {A}")
    print(f"And   B = {B}")
    print(f"The rotation index is i = {i}.")
    print("\nThis means B = A[i:] + A[:i]. Let's verify:")
    
    # "output each number in the final equation"
    rotated_A_part1 = A[i:]
    rotated_A_part2 = A[:i]
    print(f"Equation: {B} = {rotated_A_part1} + {rotated_A_part2}")
    
    # The verification is implicitly true due to the problem statement.
    # We can assert this for our own confidence.
    assert B == rotated_A_part1 + rotated_A_part2

    return i

# --- Example Usage ---
# Example 1
print("--- Example 1 ---")
A1 = [10, 20, 30, 40, 50]
B1 = [30, 40, 50, 10, 20]
find_rotation_index(A1, B1)

# Example 2 (no rotation)
print("\n--- Example 2 ---")
A2 = [1, 2, 3, 4]
B2 = [1, 2, 3, 4]
find_rotation_index(A2, B2)