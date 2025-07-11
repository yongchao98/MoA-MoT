def find_rotation_index_and_explain(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i] and explains the steps.

    This optimal O(n) algorithm works by recognizing that the first element
    of the rotated list B must correspond to the element at index 'i' in the
    original list A.

    Args:
        A: A list of n unique integers.
        B: A list which is a cyclic shift of A.
    
    Returns:
        The rotation index i.
    """
    if not A or not B:
        print("Lists cannot be empty.")
        if len(A) == len(B) == 0:
            return 0
        return -1

    if len(A) != len(B):
        print("Lists must be of the same length.")
        return -1

    # The first element of list B gives us the value A[i].
    target_value = B[0]

    # Find the index of this value in list A.
    # The list.index() method performs a linear search, which takes O(n) time.
    # This search is the performance bottleneck, and its O(n) complexity is
    # provably the best possible for an unsorted list.
    try:
        rotation_index = A.index(target_value)
    except ValueError:
        # According to the problem statement, B is always a rotation of A,
        # so this error should not occur.
        print(f"Error: Element {target_value} from B not found in A.")
        return -1

    # Print the explanation of the "equation" used to find i.
    print(f"The given lists are:")
    print(f"A = {A}")
    print(f"B = {B}\n")
    print(f"We want to find the index 'i' such that B is A rotated by 'i'.")
    print(f"The equation for the rotation is: B[0] = A[i].\n")
    
    print(f"Step 1: Get the first element of B.")
    print(f"   B[0] = {B[0]}\n")
    
    print(f"Step 2: Find the index of this element ({B[0]}) in A.")
    print(f"   The value {target_value} is located at index {rotation_index} in A, so A[{rotation_index}] = {A[rotation_index]}.\n")
    
    print(f"Step 3: Conclude the rotation index.")
    print(f"   From B[0] = A[{rotation_index}], we determine that the rotation index is i = {rotation_index}.")

    return rotation_index

# --- Example Usage ---
# Define the lists A and B
A = [15, 18, 2, 3, 6, 12]
B = [2, 3, 6, 12, 15, 18]

# Call the function to find and explain the rotation index
find_rotation_index_and_explain(A, B)