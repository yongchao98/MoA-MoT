import collections

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A: A list with n unique integers.
    B: A list which is a cyclic shift of A (B = A[i:] + A[:i]).
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: lists must be non-empty and of the same length.")
        return None

    # The core logic relies on the fact that B[0] must be A[i].
    first_element_of_B = B[0]
    print(f"The first element of list B is {first_element_of_B}.")
    print(f"We need to find the index 'i' in list A such that A[i] = {first_element_of_B}.")
    
    # Finding the index of an element in an unordered list is an O(n) operation.
    # This is the dominant part of the algorithm and has a lower bound of Omega(n).
    try:
        # A.index() performs a linear search, which is O(n).
        i = A.index(first_element_of_B)
        
        # We found the index. Now we print the result based on our core equation.
        print(f"Found the element {first_element_of_B} in list A at index {i}.")
        print("\nThe defining equation is A[i] = B[0].")
        print("Substituting the values, we get:")
        print(f"A[{i}] = {A[i]}")
        print(f"B[0] = {B[0]}")
        print(f"\nThus, the rotation index is i = {i}.")

    except ValueError:
        # This case is ruled out by the problem description but is good practice to handle.
        print(f"Error: Element {first_element_of_B} not found in list A.")
        return None

# Example usage:
A = [8, 9, 10, 11, 15, 2, 4, 6]
# Let's rotate A by i=5
# A[5:] = [2, 4, 6]
# A[:5] = [8, 9, 10, 11, 15]
# B = [2, 4, 6, 8, 9, 10, 11, 15]
B = A[5:] + A[:5]

find_rotation_index(A, B)
