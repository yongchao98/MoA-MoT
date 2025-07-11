def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A: A list with n unique integers.
    B: A list which is a cyclic shift of A (B = A[i:] + A[:i]).
    """
    n = len(A)
    # Handle empty lists as a base case.
    if n == 0:
        print("Input lists are empty, rotation index is considered 0.")
        return

    # The problem guarantees that B is a rotation of A.
    # Therefore, B[0] must be equal to A[i].
    # Since all elements in A are unique, finding the index of B[0]
    # in A will give us the unique rotation index i.
    first_element_of_B = B[0]

    # The .index() method in Python performs a linear search, which has a
    # time complexity of O(n). This aligns with our finding that O(n) is
    # the best possible complexity.
    try:
        i = A.index(first_element_of_B)
    except ValueError:
        # This case should not be reached given the problem's constraints.
        print(f"Error: B is not a valid rotation of A, as {first_element_of_B} is not in A.")
        return

    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"The rotation index is: i = {i}")

    # Output the final equation B = A[i:] + A[:i] showing each number
    part1 = A[i:]
    part2 = A[:i]
    print("\nThe verification equation is: B = A[i:] + A[:i]")
    
    # We build the string representation of the equation
    # The str() function on a list correctly formats it with brackets and commas.
    equation = f"{str(B)} = {str(part1)} + {str(part2)}"
    print(equation)


if __name__ == '__main__':
    # Example usage:
    A = [8, 9, 10, 1, 2, 3, 4, 5, 6, 7]
    # Let's create B by rotating A by i=3
    i = 3
    B = A[i:] + A[:i] # B will be [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    find_rotation_index(A, B)