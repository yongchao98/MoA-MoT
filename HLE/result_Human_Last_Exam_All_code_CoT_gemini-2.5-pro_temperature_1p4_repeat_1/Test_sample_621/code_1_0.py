def find_rotation_index():
    """
    Finds the rotation index i for two lists A and B, where B is a rotation of A.
    The time complexity is O(n) due to the search for B[0] in A.
    This is optimal as any algorithm must at least inspect all elements in the worst case.
    """
    # Let A be a list with n unique integers
    A = [3, 8, 1, 9, 4, 7, 2, 6, 5]
    n = len(A)
    
    # Let i be the rotation index. For this example, let's choose i = 5
    i_actual = 5
    
    # B is defined as A rotated by i, i.e., B = A[i:] + A[:i]
    B = A[i_actual:] + A[:i_actual]

    print(f"Given A: {A}")
    print(f"Given B: {B}")
    print("-" * 20)
    
    # The relationship B[0] == A[i] must hold.
    # We can find i by finding the index of B[0] in A.
    # List A is unsorted, so we must perform a linear search.
    # In Python, list.index() does this with O(n) complexity.
    
    if not B:
        print("List B is empty, cannot determine rotation index.")
        return

    target_element = B[0]
    
    try:
        # The search operation is the dominant part, taking O(n) time.
        found_i = A.index(target_element)
        
        print(f"The first element of B is {B[0]}.")
        print(f"Searching for {B[0]} in A gives the index {found_i}.")
        print(f"Therefore, the rotation index is i = {found_i}.")
        print("\nThe final equation is B = A[i:] + A[:i], which is:")

        # Outputting each number in the final equation by printing the lists.
        # This confirms that our found_i correctly reconstructs B from A.
        reconstructed_A_part1 = A[found_i:]
        reconstructed_A_part2 = A[:found_i]
        
        print(f"{B} = {reconstructed_A_part1} + {reconstructed_A_part2}")

    except ValueError:
        # This part should not be reached given the problem's constraints.
        print(f"Error: Element {target_element} not found in list A.")


find_rotation_index()
<<<A>>>