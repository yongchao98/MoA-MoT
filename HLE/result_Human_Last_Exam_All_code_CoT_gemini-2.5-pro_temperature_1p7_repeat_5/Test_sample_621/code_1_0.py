def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A and B are lists of unique integers.
    B is guaranteed to be a rotation of A.
    
    The time complexity is O(n) due to the call to A.index(), which performs a linear search.
    An alternative O(n) implementation would be to build a hash map of A's elements
    to their indices, which also takes O(n) time for setup.
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: lists must be non-empty and of the same size.")
        return

    # The first element of B must correspond to A[i].
    # Since all elements are unique, we just need to find the index
    # of B[0] in A.
    try:
        first_element_of_B = B[0]
        i = A.index(first_element_of_B)

        print(f"Given A = {A}")
        print(f"Given B = {B}")
        print(f"The rotation index is i = {i}")
        print("\nThis means that B is formed by taking the slice of A from index i to the end,")
        print("and concatenating it with the slice of A from the beginning up to index i.")
        print("\nEquation: B = A[i:] + A[:i]")
        
        # Demonstrating the final equation with the numbers
        # Python's list representation is used to show the numbers
        print("Final Equation with numbers:")
        a_part1 = A[i:]
        a_part2 = A[:i]
        
        # Using f-strings to print the lists which contain the numbers.
        # This fulfills the requirement to output each number in the final equation.
        print(f"{B} = {a_part1} + {a_part2}")
        
    except ValueError:
        print(f"Error: The first element of B ({B[0]}) was not found in A.")
        print("This violates the problem assumption that B is a rotation of A.")


# Example usage:
A = [8, 9, 10, 1, 2, 3, 4, 5, 6, 7]
# Let's rotate A by i=3 to create B
# A[3:] = [1, 2, 3, 4, 5, 6, 7]
# A[:3] = [8, 9, 10]
# B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
B = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

find_rotation_index(A, B)
