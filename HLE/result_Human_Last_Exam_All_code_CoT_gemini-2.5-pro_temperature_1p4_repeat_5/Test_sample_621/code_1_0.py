def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    The algorithm has a time complexity of O(n) because it requires searching for
    an element in an unsorted list. This is proven to be the optimal
    complexity because any algorithm must, in the worst case, scan the entirety
    of list A to locate the position of the first element of B.

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The integer rotation index i.
    """
    if not A or not B:
        print("Input lists cannot be empty.")
        return None

    # The first element of the rotated list B must be A[i].
    target_element = B[0]
    n = len(A)
    rotation_index = -1

    # Search for the target element in A. This is the O(n) step.
    for i in range(n):
        if A[i] == target_element:
            rotation_index = i
            break
    
    if rotation_index == -1:
        # This case should not be reached given the problem's constraints.
        print("B does not appear to be a rotation of A.")
        return None

    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"\nStep 1: The first element of B is {target_element}.")
    print(f"Step 2: Searching for {target_element} in A.")
    print(f"Step 3: Found {target_element} at index {rotation_index} in A.")
    print(f"\nConclusion: The rotation index is i = {rotation_index}.")

    # Verification and printing the "final equation"
    # This part is for demonstration and doesn't affect the complexity.
    reconstructed_b = A[rotation_index:] + A[:rotation_index]
    
    print("\nFinal Equation:")
    print(f"B = A[i:] + A[:i]")
    
    # "output each number in the final equation!"
    # We interpret this as printing the lists and the index i.
    b_str = ', '.join(map(str, B))
    reconstructed_b_str = ', '.join(map(str, reconstructed_b))

    print(f"[{b_str}] = A[{rotation_index}:] + A[:{rotation_index}]")
    print(f"Which evaluates to:")
    print(f"[{b_str}] = [{reconstructed_b_str}]")


# --- Example Usage ---
# Define sample lists A and B
list_A = [15, 28, 5, 42, 33, 19, 50, 2]
# Let's create B by rotating A with i=3
# A[3:] = [42, 33, 19, 50, 2]
# A[:3] = [15, 28, 5]
# B = [42, 33, 19, 50, 2, 15, 28, 5]
list_B = [42, 33, 19, 50, 2, 15, 28, 5]

find_rotation_index(list_A, list_B)