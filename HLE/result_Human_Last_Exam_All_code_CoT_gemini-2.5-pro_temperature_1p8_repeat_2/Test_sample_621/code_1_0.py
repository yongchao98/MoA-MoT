def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    Args:
        A: A list with n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i, or -1 if B is empty or its first element is not in A.
    """
    if not B:
        print("Error: List B is empty.")
        return -1

    # The first element of the rotated list B must be A[i].
    target_element = B[0]

    # Find the index of the target element in list A.
    # In an unsorted list, this search operation has a worst-case time complexity of O(n).
    try:
        # Using list.index() is a clean way to perform the search.
        found_index = A.index(target_element)
        return found_index
    except ValueError:
        # This case should not occur given the problem statement.
        print(f"Error: Element {target_element} not found in list A.")
        return -1

def solve():
    """
    Sets up an example and calls the function to solve the problem.
    """
    A = [15, 25, 35, 45, 55, 65, 75]
    # Let's choose a rotation index for our example
    i_actual = 4
    B = A[i_actual:] + A[:i_actual]

    print(f"Given list A: {A}")
    print(f"Given list B (a rotation of A): {B}")
    print("-" * 20)

    found_i = find_rotation_index(A, B)

    if found_i is not None and found_i != -1:
        print(f"The algorithm finds that the first element of B is {B[0]}.")
        print(f"Searching for {B[0]} in A, we find it at index {found_i}.")
        print(f"Therefore, the rotation index is i = {found_i}.")
        
        # As requested, output the final equation showing all numbers
        print("\nFinal Equation Verification:")
        print(f"B                         = {B}")
        # Construct the string for the right-hand side of the equation
        rhs_expr = f"A[{found_i}:] + A[:{found_i}]"
        rhs_val = A[found_i:] + A[:found_i]
        print(f"{rhs_expr:25} = {A[found_i:]} + {A[:found_i]}")
        print(f"{'':25} = {rhs_val}")

        if B == rhs_val:
            print("\nThe equation holds true.")
        else:
            print("\nError: The found index is incorrect.")


solve()