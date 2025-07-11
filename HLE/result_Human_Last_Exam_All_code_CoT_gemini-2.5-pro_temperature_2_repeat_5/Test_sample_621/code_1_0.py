def solve():
    """
    Finds the rotation index `i` for two lists A and B, where B is a rotation of A.
    The algorithm runs in O(n) time, which is optimal.
    """
    # Example lists. You can replace these with any other valid inputs.
    A = [3, 8, 1, 9, 4, 2, 7, 5]
    true_i = 4
    B = A[true_i:] + A[:true_i]  # Creates the rotated list: [4, 2, 7, 5, 3, 8, 1, 9]

    print(f"Given list A: {A}")
    print(f"Given list B: {B}")
    print("-" * 25)

    if len(A) != len(B):
        print("Error: Lists must have the same length.")
        return
    if not A:
        print("Lists are empty. The rotation index is 0.")
        i = 0
        print(f"\nThe final equation is:")
        print(f"B = A[{i}:] + A[:{i}]")
        print(f"[] = [] + []")
        return

    # The key insight is that B[0] must equal A[i]. Since all elements in A are unique,
    # we can find 'i' by simply locating the first element of B within A.
    target_element = B[0]

    try:
        # Searching for the element in list A. This operation has a time complexity of O(n).
        i = A.index(target_element)
    except ValueError:
        # This case would not happen if B is guaranteed to be a rotation of A.
        print(f"Error: Element {target_element} not found in A. B cannot be a rotation of A.")
        return
    
    print(f"The first element of B is {target_element}.")
    print(f"This element is found at index {i} in list A.")
    print(f"Therefore, the rotation index is {i}.")

    # As requested, output the equation with the actual numbers.
    print(f"\nThe final equation is:")
    
    # Represent the lists and the operation
    a_part1_str = str(A[i:])
    a_part2_str = str(A[:i])
    b_str = str(B)
    
    print(f"B = A[{i}:] + A[:{i}]")
    print(f"{b_str} = {a_part1_str} + {a_part2_str}")

# Execute the solver function
solve()