def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is A rotated by i.
    Time complexity: O(n) due to the A.index() call.
    Space complexity: O(1) (or O(n) for slices if you count output, but core logic is O(1)).
    """
    if not A or len(A) != len(B):
        print("Invalid input. Lists must be non-empty and of the same length.")
        return

    n = len(A)
    # The first element of B must exist in A. Since all elements are unique,
    # its index in A determines the rotation amount.
    # The .index() method has a time complexity of O(n).
    try:
        i = A.index(B[0])
    except ValueError:
        # This case is excluded by the problem statement, which guarantees
        # that B is a rotation of A.
        print("B is not a rotation of A.")
        return

    # To satisfy the "output each number in the final equation" request,
    # we demonstrate the result with the actual list slices.
    print(f"Given A = {A}")
    print(f"And B = {B}")
    print(f"The rotation index is i = {i}.")
    print("\nThis means B = A[i:] + A[:i]. Let's verify:")
    
    rotated_part1 = A[i:]
    rotated_part2 = A[:i]
    
    print(f"A[i:] (i={i}) is A[{i}:] = {rotated_part1}")
    print(f"A[:i] (i={i}) is A[:{i}] = {rotated_part2}")
    
    reconstructed_B = rotated_part1 + rotated_part2
    print(f"A[{i}:] + A[:{i}] = {reconstructed_B}")
    
    if reconstructed_B == B:
        print("The reconstructed list matches B.")
    else:
        # This should not happen given the problem constraints.
        print("The reconstructed list does NOT match B.")

# Example Usage:
# Let A be a list of n unique integers.
A = [10, 20, 30, 40, 50, 60]
# Let B be a rotation of A, for example with i=2
# B = A[2:] + A[:2] = [30, 40, 50, 60] + [10, 20]
B = [30, 40, 50, 60, 10, 20]

find_rotation_index(A, B)
