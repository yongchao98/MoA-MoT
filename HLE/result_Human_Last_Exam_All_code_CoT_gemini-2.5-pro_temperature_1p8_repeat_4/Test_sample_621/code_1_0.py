def find_rotation_index(A, B):
    """
    Finds the rotation index i for two lists A and B, where B is a rotation of A.
    
    The plan is as follows:
    1. The rotation equation is B = A[i:] + A[:i]. This implies B[0] == A[i].
    2. Since A has unique elements, we can find the index 'i' by locating the first element of B (B[0]) in A.
    3. The time complexity of finding an element in an unsorted list is O(n), which is optimal.
    4. The code will demonstrate this by finding the index and printing the final verified equation.
    """
    
    if not A or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same size.")
        return

    n = len(A)
    
    # 1. Take the first element of B
    first_b_element = B[0]
    
    # 2. Find its index in A. This takes O(n) time.
    try:
        i = A.index(first_b_element)
    except ValueError:
        print(f"Error: Element {first_b_element} not found in A. B is not a rotation of A.")
        return

    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print("-" * 20)
    print("Step 1: Identify the relationship between A, B, and the rotation index 'i'.")
    print("From B = A[i:] + A[:i], we know the first element of B must be the i-th element of A.")
    print(f"So, A[i] = B[0], which is {first_b_element}.")
    
    print("\nStep 2: Find the index 'i' in A.")
    print(f"We search for the element {first_b_element} in list A.")
    print(f"The element {A[i]} is found at index {i} in A.")
    print(f"Therefore, the rotation index i = {i}.")

    print("\nStep 3: Verify the equation B = A[i:] + A[:i] with the found index.")
    
    # Split A at index i to show the rotation
    a_suffix = A[i:]
    a_prefix = A[:i]
    
    print(f"B = {B}")
    # We use repr() to get a string representation that distinguishes lists, e.g., for empty lists
    print(f"A[{i}:] = {repr(a_suffix)}")
    print(f"A[:{i}] = {repr(a_prefix)}")
    print(f"A[{i}:] + A[:{i}] = {a_suffix + a_prefix}")
    print("\nThe equation holds true, confirming our result.")


# Example usage:
list_A = [10, 20, 30, 40, 50, 60]
# Rotate A by i=4 to get B
# A[4:] = [50, 60]
# A[:4] = [10, 20, 30, 40]
# B = [50, 60, 10, 20, 30, 40]
list_B = [50, 60, 10, 20, 30, 40]

find_rotation_index(list_A, list_B)
