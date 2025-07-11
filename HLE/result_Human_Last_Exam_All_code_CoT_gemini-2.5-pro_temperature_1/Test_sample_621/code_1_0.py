def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i.
    """
    if not A or not B or len(A) != len(B):
        print("Error: Lists must be non-empty and of the same length.")
        return -1

    if len(A) == 0:
        return 0
        
    n = len(A)

    # The first element of B must be A[i].
    # We can find i by finding the index of B[0] in A.
    # This step has a time complexity of O(n) because list.index()
    # performs a linear scan in the worst case.
    pivot = B[0]
    try:
        i = A.index(pivot)
    except ValueError:
        print(f"Error: {pivot} not found in A. B is not a cyclic shift of A.")
        return -1

    # Now we have the index i. We can print it and verify the result.
    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"The found rotation index is: i = {i}")

    # To demonstrate the result, we show the equation B = A[i:] + A[:i]
    # with the actual numbers.
    rotated_A = A[i:] + A[:i]
    print(f"\nVerification:")
    print(f"B = {B}")
    print(f"A[i:] = {A[i:]}")
    print(f"A[:i] = {A[:i]}")
    print(f"A[i:] + A[:i] = {rotated_A}")
    
    if B == rotated_A:
        print("\nThe lists match, the index is correct.")
    else:
        print("\nVerification failed. The index is incorrect.")
        
    return i

# Example usage:
# Let A = [10, 20, 30, 40, 50] and i = 3
# Then B should be A[3:] + A[:3] = [40, 50] + [10, 20, 30] = [40, 50, 10, 20, 30]
A_list = [10, 20, 30, 40, 50]
B_list = [40, 50, 10, 20, 30]

find_rotation_index(A_list, B_list)