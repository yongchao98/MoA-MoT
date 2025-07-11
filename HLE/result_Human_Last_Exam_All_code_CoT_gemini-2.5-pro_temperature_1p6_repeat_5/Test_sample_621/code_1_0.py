def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a rotation of A.
    B = A[i:] + A[:i]
    The time complexity is O(n) because of the index search and list comparison.

    Args:
        A: A list with n unique integers.
        B: A list that is potentially a rotation of A.

    Returns:
        The integer index i if B is a rotation of A, otherwise -1.
    """
    n = len(A)
    if n != len(B):
        return -1
    if n == 0:
        return 0

    # The first element of B must be somewhere in A.
    # Find the index of B[0] in A. This is O(n).
    try:
        first_b_in_a_index = A.index(B[0])
    except ValueError:
        # If the first element of B is not in A, it's not a rotation.
        return -1
    
    # We found a candidate index i. Now, verify if the rotation matches B.
    # This list slicing and concatenation creates a new list, taking O(n) time and space.
    # The subsequent comparison also takes O(n).
    if A[first_b_in_a_index:] + A[:first_b_in_a_index] == B:
        return first_b_in_a_index
    else:
        return -1

# Example usage:
# A is a list of n unique integers.
# B is A rotated by i=17
A = [20, 27, 3, 11, 42, 18, 7, 33, 1, 29, 38, 14, 49, 22, 10, 31, 45, 5, 25, 40]
B = [5, 25, 40, 20, 27, 3, 11, 42, 18, 7, 33, 1, 29, 38, 14, 49, 22, 10, 31, 45]

# Find the rotation index 'i'
rotation_index = find_rotation_index(A, B)

# Per instructions, print the numbers involved. 
# Here, the main "number" in the final result is the index 'i'.
if rotation_index != -1:
    print(f"Given A: {A}")
    print(f"Given B: {B}")
    print(f"Found rotation index i = {rotation_index}")
    # Verification of the "equation" B = A[i:] + A[:i]
    print(f"A rotated by {rotation_index} is: {A[rotation_index:] + A[:rotation_index]}")
else:
    print("B is not a rotation of A.")