def find_rotation_index(A, B):
    """
    Given two lists A and B, where B is a cyclic shift of A (A[i:] + A[:i]),
    this function finds the rotation index i.

    The problem states that A contains n unique integers and is not necessarily ordered.
    B is guaranteed to be a rotation of A.

    The plan is as follows:
    1. The relationship B = A[i:] + A[:i] implies that B[0] must equal A[i].
    2. Therefore, the task is to find the index of the element B[0] inside list A.
    3. Since A is an unsorted list, searching for an element requires a linear scan,
       which has a time complexity of O(n).
    4. This O(n) approach is optimal because the lower bound for searching in an
       unsorted array is Ω(n).

    Args:
        A (list): A list of n unique integers.
        B (list): A list that is a rotation of A.

    Returns:
        int: The rotation index i.
    """
    if not A or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same length.")
        return -1

    # Get the target element, which is the first element of the rotated list B.
    target = B[0]

    try:
        # Use list.index() to find the first occurrence of the target in A.
        # This operation is O(n) in the worst case.
        i = A.index(target)
        
        print(f"Given A: {A}")
        print(f"Given B: {B}")
        print(f"The first element of B is {target}.")
        print(f"This element is located at index {i} in list A.")
        # As per the prompt, we output the number in the 'final equation'.
        # We interpret this as printing the final resulting index 'i'.
        print(f"The rotation index 'i' is: {i}")

        return i
    except ValueError:
        # This block will not be reached if B is guaranteed to be a rotation of A.
        print(f"Error: Element {target} from B was not found in A.")
        return -1

# Example usage of the function
# Let A be a list of unique, unsorted integers
A_list = [23, 51, 19, 8, 42, 65]
# Let B be A rotated by i = 2
# A[2:] = [19, 8, 42, 65]
# A[:2] = [23, 51]
# B = [19, 8, 42, 65, 23, 51]
B_list = [19, 8, 42, 65, 23, 51]

# Call the function to find the rotation index
find_rotation_index(A_list, B_list)