def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i.
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: lists must be non-empty and of the same length.")
        return -1

    # The first element of B must be A[i].
    element_to_find = B[0]

    # Find the index of this element in A.
    # This is a linear search, which has O(n) time complexity.
    try:
        i = A.index(element_to_find)
    except ValueError:
        print(f"Error: Element {element_to_find} from B was not found in A.")
        return -1
        
    print(f"Given list A: {A}")
    print(f"Given list B: {B}")
    print("-" * 20)
    print(f"The first element of B is: {B[0]}")
    # The prompt asks to "output each number in the final equation".
    # We can interpret this as showing the logic.
    # The "equation" is finding i such that A[i] == B[0].
    print(f"We need to find the index 'i' in A where A[i] == {B[0]}.")
    print(f"The element {B[0]} is found at index {i} in A.")
    print(f"Therefore, the rotation index is: {i}")
    
    # Verification step
    reconstructed_B = A[i:] + A[:i]
    print(f"Verification: A[{i}:] + A[:{i}] = {reconstructed_B}")
    if reconstructed_B == B:
        print("Verification successful.")
    else:
        print("Verification failed.")
        
    return i

# Example Usage:
# Let's define A and create B by rotating A with i=6
A = [15, 23, 8, 42, 16, 4, 99, 54, 77, 2]
i_original = 6
B = A[i_original:] + A[:i_original]

# Now, let's find i using our algorithm
find_rotation_index(A, B)
