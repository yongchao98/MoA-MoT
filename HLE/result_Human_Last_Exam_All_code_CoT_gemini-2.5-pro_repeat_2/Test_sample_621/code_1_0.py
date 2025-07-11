def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    The list A is assumed to contain unique integers.

    This implementation has a time complexity of O(n) and space complexity of O(n).

    Args:
        A: A list of n unique integers (not necessarily ordered).
        B: A list which is a cyclic shift of A.

    Returns:
        The rotation index i, or -1 if inputs are invalid.
    """
    if not A or not B or len(A) != len(B):
        print("Error: Input lists must be non-empty and of the same size.")
        return -1
    
    n = len(A)
    
    # For an efficient O(1) average-time lookup, we create a map from
    # each value in A to its index. This step takes O(n) time.
    value_to_index_map = {value: i for i, value in enumerate(A)}
    
    # The first element of B must be A[i]. We get this target value.
    # This is an O(1) operation.
    target_value = B[0]
    
    # We find the index of the target value in A using our map.
    # This lookup is O(1) on average.
    if target_value not in value_to_index_map:
        # This case should not be reached if B is truly a rotation of A.
        print(f"Error: B cannot be a rotation of A, because the element {target_value} is not in A.")
        return -1
        
    rotation_index = value_to_index_map[target_value]
    
    # The problem guarantees that B is a rotation of A and A has unique elements,
    # so this index is guaranteed to be the correct rotation index `i`.
    
    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"\nThe rotation index 'i' is: {rotation_index}")
    
    # Displaying the final equation with all its numbers as requested
    print("\nThe equation B = A[i:] + A[:i] with the found index is:")
    print(f"{B} = {A[rotation_index:]} + {A[:rotation_index]}")


# --- Example Usage ---
# Let's define our lists A and B.
# Here, B is A rotated by i=3.
A = [15, 28, 32, 49, 51, 67]
B = [49, 51, 67, 15, 28, 32]

# Run the algorithm and print the result.
find_rotation_index(A, B)