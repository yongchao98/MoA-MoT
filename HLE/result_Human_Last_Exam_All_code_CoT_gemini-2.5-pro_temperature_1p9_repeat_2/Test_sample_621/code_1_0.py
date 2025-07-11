def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This function relies on the fact that A contains unique elements and B is a 
    guaranteed rotation of A. This means that the first element of B, B[0], 
    must correspond to the element A[i]. By finding the index of B[0] in A, 
    we can determine i.

    The time complexity is O(n) because finding the index of an element in a list
    requires, in the worst case, scanning the entire list. As explained in the analysis,
    this is the optimal time complexity for this problem.
    """
    # Check for empty lists as an edge case.
    if not A and not B:
        print("Input lists are empty. Rotation index is 0.")
        print("[] = [] + []")
        return 0
    if not A or not B or len(A) != len(B):
        print("Invalid input: lists must be non-empty and have the same size.")
        return -1

    # The first element of the rotated list B must be A[i].
    target_element = B[0]

    # Find the index of this element in A.
    # list.index() in Python has a time complexity of O(n).
    try:
        i = A.index(target_element)
    except ValueError:
        # This case should not be reached if B is guaranteed to be a rotation of A.
        print(f"Error: Element {target_element} not found in A.")
        return -1
    
    print(f"Given A: {A}")
    print(f"Given B: {B}")
    print(f"The rotation index 'i' is: {i}")
    
    # Reconstruct the rotation equation and print it.
    rotated_part1 = A[i:]
    rotated_part2 = A[:i]
    
    # To prevent large lists from printing on a single line and wrapping messily,
    # we convert lists to strings for printing.
    str_B = ", ".join(map(str, B))
    str_p1 = ", ".join(map(str, rotated_part1))
    str_p2 = ", ".join(map(str, rotated_part2))
    
    print(f"The equation B = A[i:] + A[:i] is:")
    print(f"[{str_B}] = [{str_p1}] + [{str_p2}]")

    return i

# --- Example Usage ---
# Let A be a list with n unique integers.
A = [8, 9, 10, 1, 2, 3, 4, 5, 6, 7]
# Let B be a rotation of A by i=3.
# A[3:] = [1, 2, 3, 4, 5, 6, 7]
# A[:3] = [8, 9, 10]
# B = A[3:] + A[:3]
B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Run the algorithm to find the index 'i'.
find_rotation_index(A, B)
