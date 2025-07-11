def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    Args:
        A: A list of n unique integers.
        B: A list which is a rotation of A.

    Returns:
        The rotation index i, or None if B is not a rotation of A.
    """
    n = len(A)
    if n != len(B):
        print("B is not a rotation of A because lists are of different lengths.")
        return
    if n == 0:
        print("The rotation index i is: 0") # By convention for empty lists
        return

    # Take the first element of B to find the potential starting point in A.
    # Uniqueness in A guarantees at most one candidate for i.
    first_b_element = B[0]

    try:
        # list.index() performs a linear scan, taking O(n) time.
        candidate_i = A.index(first_b_element)
    except ValueError:
        # If the first element of B is not in A, it can't be a rotation.
        print("B is not a rotation of A.")
        return

    # Verify that B is indeed A rotated by candidate_i.
    # This loop runs n times, taking O(n) time.
    # It compares elements without creating a new rotated list in memory.
    for k in range(n):
        if B[k] != A[(candidate_i + k) % n]:
            print("B is not a rotation of A.")
            return

    # If the loop completes, the verification is successful.
    # The prompt asks to output the number, which we interpret as printing the index i.
    print(f"The rotation index i is: {candidate_i}")


# Example usage:
A = [10, 20, 30, 40, 50, 60]
# B is A rotated by i=2
B = [30, 40, 50, 60, 10, 20]

find_rotation_index(A, B)
