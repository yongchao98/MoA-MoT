def find_rotation_index(A, B):
    """
    Given two lists A and B, where B is a cyclic shift of A, this function
    finds the rotation index i such that B = A[i:] + A[:i].
    The lists are assumed to contain unique integers.
    The algorithm has a time complexity of O(n).
    """
    if not A or not B or len(A) != len(B):
        print("Error: Input lists must be non-empty and have the same size.")
        return

    n = len(A)
    # The first element of the rotated list B, which is B[0], corresponds to A[i].
    # Because all elements in A are unique, finding the index of B[0] in A
    # will uniquely identify the rotation index i.
    element_to_find = B[0]

    try:
        # The list.index() method performs a linear scan, which has O(n) complexity.
        # This is the dominant part of the algorithm.
        i = A.index(element_to_find)
    except ValueError:
        # This case should not happen if B is guaranteed to be a rotation of A.
        print(f"Error: The element {element_to_find} (from B[0]) was not found in A.")
        return

    print(f"Given list A: {A}")
    print(f"Given list B: {B}")
    print(f"The rotation index 'i' is the index of B[0] (which is {B[0]}) in A.")
    print(f"We found that A[{i}] = {A[i]}.")
    print(f"Therefore, the rotation index is i = {i}.")
    print("-" * 20)
    
    # As requested, printing the final equation with all numbers.
    print("Verification of the equation B = A[i:] + A[:i]:")
    rotated_A_part1 = A[i:]
    rotated_A_part2 = A[:i]
    
    print(f"B = {B}")
    print(f"A[{i}:] = {rotated_A_part1}")
    print(f"A[:{i}] = {rotated_A_part2}")
    print(f"A[i:] + A[:i] = {rotated_A_part1 + rotated_A_part2}")
    # Final confirmation
    is_correct = (B == rotated_A_part1 + rotated_A_part2)
    print(f"\nDoes B equal A[i:] + A[:i]? {is_correct}")

# --- Example Usage ---
# Define the original list A.
A = [1, 2, 3, 4, 5, 6, 7]
# Define the list B, which is A rotated by i=4.
# A[4:] = [5, 6, 7]
# A[:4] = [1, 2, 3, 4]
# B = [5, 6, 7, 1, 2, 3, 4]
B = [5, 6, 7, 1, 2, 3, 4]

find_rotation_index(A, B)