def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This function assumes:
    - A and B are lists of the same length n.
    - A contains unique integers.
    - B is a cyclic shift of A.

    Complexity: O(n) because of the .index() call and list slicing/comparison.
    This is optimal as any algorithm must verify all n elements in the worst case.
    """
    if not A or not B or len(A) != len(B):
        print("Input lists must be non-empty and have the same length.")
        return

    n = len(A)
    if n == 0:
        print("Lists are empty, rotation index is 0.")
        return

    # Find the first element of B in A. This gives the rotation index.
    # The .index() method takes O(n) time.
    first_element_of_B = B[0]
    try:
        i = A.index(first_element_of_B)
    except ValueError:
        print(f"Element {first_element_of_B} not found in A. B is not a rotation of A.")
        return

    # Verify that the rotation is correct. This step is implicitly O(n)
    # due to slicing and comparison, and necessary for correctness proof,
    # though with the problem's guarantees, it will always pass.
    # If A_rotated == B, then we have found our index.
    A_rotated = A[i:] + A[:i]
    if A_rotated == B:
        print(f"Given A = {A}")
        print(f"Given B = {B}")
        print(f"The rotation index is i = {i}.")
        print("\nFinal Equation:")
        # The prompt requires printing each number in the final equation.
        # This format shows the lists involved.
        print(f"{B} = {A[i:]} + {A[:i]}")
    else:
        # This case won't be reached if B is guaranteed to be a rotation of A.
        print("Error: B is not a simple cyclic shift of A.")


# Example Usage:
list_A = [8, 9, 10, 1, 2, 3, 4, 5, 6, 7]
list_B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
find_rotation_index(list_A, list_B)