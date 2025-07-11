def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This function has a time complexity of O(n) and is optimal.

    Args:
        A: The original list of n unique integers.
        B: The rotated list.

    Returns:
        The rotation index i.
    """
    if not A or not B or len(A) != len(B):
        # Basic validation for empty lists or lists of different sizes
        return -1

    n = len(A)
    # The list A contains unique integers. The first element of B, B[0],
    # must be present somewhere in A. Finding its index gives us the rotation value.
    # The 'index' method performs a linear search, which is O(n).
    try:
        first_element_b = B[0]
        i = A.index(first_element_b)
    except ValueError:
        # This case happens if B[0] is not in A, which means B is not a rotation of A.
        return -1
        
    # The problem statement guarantees B is a rotation of A, so we don't need
    # to perform a full verification which would be:
    # assert B == A[i:] + A[:i]
    # This verification would take an additional O(n) time but doesn't change the
    # overall complexity.

    return i

# --- Example Usage ---

# Define the lists A and B
A = [10, 20, 30, 40, 50, 60]
B = [40, 50, 60, 10, 20, 30]

print(f"Given list A: {A}")
print(f"Given list B: {B}")

# Find the rotation index i
i = find_rotation_index(A, B)

if i != -1:
    print(f"\nThe rotation index is: {i}")
    
    # To satisfy the prompt "output each number in the final equation",
    # we print the equation showing how B is constructed from A using index i.
    rotated_part1 = A[i:]
    rotated_part2 = A[:i]
    
    print("\nThe equation B = A[i:] + A[:i] is demonstrated as:")
    print(f"{B} = {rotated_part1} + {rotated_part2}")

else:
    print("B is not a valid rotation of A.")
