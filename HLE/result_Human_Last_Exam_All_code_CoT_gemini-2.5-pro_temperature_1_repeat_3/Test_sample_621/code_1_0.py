def find_rotation_index():
    """
    Given two lists A and B, where B is a cyclic shift of A, this function
    finds the rotation index i.
    """
    # Example lists. A has n unique integers.
    # B = A[i:] + A[:i]
    A = [15, 20, 25, 30, 5, 10]
    i_actual = 4
    B = A[i_actual:] + A[:i_actual]

    print(f"Goal: Find the rotation index 'i' for the lists:")
    print(f"A = {A}")
    print(f"B = {B}")
    print("-" * 20)
    
    found_i = -1
    if not A and not B:
        found_i = 0
    elif A and B:
        # The algorithm is based on the fact that B[0] must be A[i].
        # Since A contains unique elements, we can find 'i' by locating
        # the first element of B in A.
        # The list.index() method has a time complexity of O(n).
        try:
            first_element_of_B = B[0]
            found_i = A.index(first_element_of_B)
        except ValueError:
            # This case should not be reached if B is guaranteed to be a rotation of A.
            print("Error: B is not a rotation of A.")
            return

    print(f"The rotation index found by the algorithm is: i = {found_i}")

    # Outputting the equation with each number as requested.
    # The format will be B = A_part1 + A_part2
    b_str = ", ".join(map(str, B))
    
    # Check if found_i is valid before slicing
    if found_i != -1:
        rotated_part1 = A[found_i:]
        rotated_part2 = A[:found_i]
        
        a_rotated_part1_str = ", ".join(map(str, rotated_part1))
        a_rotated_part2_str = ", ".join(map(str, rotated_part2))

        print("\nVerification using the found index:")
        # Display the equation that reconstructs B from A using i.
        if rotated_part2: # If the second part is not empty
            print(f"[{b_str}] = [{a_rotated_part1_str}] + [{a_rotated_part2_str}]")
        else: # If the rotation is 0, the second part is empty
            print(f"[{b_str}] = [{a_rotated_part1_str}]")
    
find_rotation_index()