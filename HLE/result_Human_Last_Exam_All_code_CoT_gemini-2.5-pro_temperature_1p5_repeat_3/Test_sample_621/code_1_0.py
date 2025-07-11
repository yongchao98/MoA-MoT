def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    Args:
        A (list): The original list of n unique integers.
        B (list): The rotated list.

    Returns:
        int: The rotation index i.
    """
    if not A or not B or len(A) != len(B):
        print("Input lists are invalid (empty or different sizes).")
        return -1

    # The first element of the rotated list B must be the element at index i in A.
    first_element_of_B = B[0]
    
    try:
        # Find the index of this element in the original list A.
        # This is an O(n) operation as it may require a full scan of A.
        rotation_index_i = A.index(first_element_of_B)
        
        print(f"The lists are:")
        print(f"A = {A}")
        print(f"B = {B}")
        print("-" * 20)
        print("To find the shift 'i', we use the relationship: B[0] = A[i]")
        # The final equation part
        print(f"The first element of B is B[0] = {first_element_of_B}.")
        print(f"The element {first_element_of_B} is found at index {rotation_index_i} in A.")
        print(f"Therefore, the final equation is i = {rotation_index_i}.")
        
        return rotation_index_i
    except ValueError:
        print(f"Error: The element {first_element_of_B} from B was not found in A.")
        print("This suggests B is not a valid rotation of A.")
        return -1

# Example usage:
# Let A be a list of unique integers.
A = [10, 20, 30, 40, 50, 60]
# Let i = 4.
i = 4
# B is A rotated by i.
B = A[i:] + A[:i]

find_rotation_index(A, B)