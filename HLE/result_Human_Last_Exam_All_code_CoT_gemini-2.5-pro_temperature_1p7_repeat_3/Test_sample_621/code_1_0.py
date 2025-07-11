def find_rotation_index(A, B):
    """
    Given a list of unique integers A, and a list B which is a cyclic shift of A,
    this function finds the rotation index i such that B = A[i:] + A[:i].
    The time complexity of this algorithm is O(n), which is optimal.

    Args:
        A: A list of n unique integers.
        B: A list which is a rotation of A.

    Returns:
        The rotation index i.
    """
    print(f"Goal: Find index 'i' for A={A} and B={B}")
    print("-" * 30)

    # Basic validation
    if not A or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same size.")
        return -1
    
    n = len(A)
    # Step 1: The first element of B must be the i-th element of A.
    first_element_b = B[0]
    print(f"The relationship is B[0] = A[i].")
    print(f"We take the first element of B, which is B[0] = {first_element_b}.")

    # Step 2: Find the index of this element in A.
    # For an unsorted list, this is an O(n) operation.
    try:
        # list.index() performs a linear search.
        i = A.index(first_element_b)
        
        # Step 3: Because A has unique elements and B is a guaranteed rotation,
        # this index 'i' is the one we're looking for.
        print(f"We search for {first_element_b} in A.")
        print(f"It is found at index {i}, since A[{i}] = {A[i]}.")
        print(f"\nThe equation B[0] = A[i] becomes:")
        print(f"{B[0]} = {A[i]}")
        print(f"\nTherefore, the rotation index is i = {i}.")
        
        return i
    except ValueError:
        # This case should not be reached given the problem's constraints.
        print(f"Error: {first_element_b} was not found in A.")
        return -1

# Example usage:
# Define the lists A and B. B is A rotated by 3 positions.
A_list = [10, 20, 30, 40, 50, 60]
true_i = 3
B_list = A_list[true_i:] + A_list[:true_i]  # Creates [40, 50, 60, 10, 20, 30]

# Run the algorithm and print the result.
find_rotation_index(A_list, B_list)