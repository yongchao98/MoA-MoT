def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A and B are lists of n unique integers. B = A[i:] + A[:i].
    
    The time complexity is O(n) because finding the index of an element
    in an unsorted list requires a linear scan in the worst case.
    This is proven to be the optimal complexity.

    Args:
        A (list): The original list of unique integers.
        B (list): The rotated list.

    Returns:
        int: The rotation index i, or -1 if inputs are invalid.
    """
    # An empty list cannot be rotated.
    if not A or not B or len(A) != len(B):
        return -1
    
    # The first element of the rotated list B must be A[i].
    # Since all elements in A are unique, there is only one such i.
    target_element = B[0]
    
    try:
        # The list.index() method performs a linear search, which is O(n).
        i = A.index(target_element)
        return i
    except ValueError:
        # This block will not be reached if B is guaranteed to be a rotation of A.
        # It's included for robustness.
        print(f"Error: The first element of B ({target_element}) was not found in A.")
        return -1

# --- Example Usage ---
# Define the original list A and the rotation index i
A = [45, 62, 12, 35, 71, 28, 53]
i_actual = 4

# Construct the rotated list B based on A and i
B = A[i_actual:] + A[:i_actual]

# Use the function to find the rotation index
found_i = find_rotation_index(A, B)

# Print the results
print(f"Original list A: {A}")
print(f"Rotated list B:  {B}")
if found_i != -1:
    print(f"The calculated rotation index is: {found_i}")
    # This line demonstrates the final equation using the found index
    print(f"B = A[{found_i}:] + A[:{found_i}]")
else:
    print("Could not determine the rotation index.")
