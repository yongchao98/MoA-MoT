def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    Args:
        A: A list with n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i.
    """
    if not A or not B or len(A) != len(B):
        # Basic validation for empty or mismatched length lists
        return -1

    # The first element of B must be the i-th element of A.
    target_value = B[0]

    # Find the index of this target value in A.
    # The .index() method performs a linear search, which is O(n).
    try:
        i = A.index(target_value)
        return i
    except ValueError:
        # This case should not happen based on the problem description,
        # but it's good practice to handle it.
        return -1

# --- Example Usage ---
# Define the original list A
A = [15, 25, 35, 45, 55, 65, 75]
n = len(A)
# Let's choose a rotation index i, for example i = 4
i_actual = 4
# Create the rotated list B = A[4:] + A[:4]
B = A[i_actual:] + A[:i_actual]

print(f"Original list A: {A}")
print(f"Rotated list B: {B}")

# Find the rotation index using the algorithm
found_i = find_rotation_index(A, B)

# Print the final result
if found_i != -1:
    print(f"\nThe algorithm needs to solve B[0] == A[i].")
    print(f"Here, B[0] = {B[0]}.")
    print(f"The value {B[0]} is found at index {found_i} in list A.")
    print(f"Therefore, the rotation index i is {found_i}.")
else:
    print("Could not find the rotation index. The lists might not be rotations of each other.")
