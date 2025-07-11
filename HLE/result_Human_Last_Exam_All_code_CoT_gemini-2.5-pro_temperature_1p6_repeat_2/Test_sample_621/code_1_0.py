import collections

def find_rotation_index(A, B):
    """
    Given a list A of unique integers and a list B, which is a rotation of A,
    this function finds the rotation index i such that B = A[i:] + A[:i].

    The time complexity of this solution is O(n), which is optimal.
    
    Args:
      A: A list of n unique integers.
      B: A list which is a known rotation of A.

    Returns:
      None. It prints the rotation index i.
    """
    print(f"Original list A: {A}")
    print(f"Rotated list B:  {B}\n")
    
    if not A or not B:
        print("One or both lists are empty.")
        if not A and not B:
            print("The rotation index for empty lists can be considered 0.")
        return

    n = len(A)
    # The element at the start of the rotated list B is B[0].
    # This element was at index `i` in the original list A.
    element_to_find = B[0]

    # To find the index `i`, we can search for `element_to_find` in `A`.
    # Using A.index() performs a linear search, which is O(n).
    try:
        i = A.index(element_to_find)
    except ValueError:
        # This part should not be reached given the problem's constraints.
        print("Error: B does not seem to be a rotation of A.")
        return
        
    print(f"The goal is to find 'i' in the equation: B = A[i:] + A[:i]")
    print(f"1. The first element of B is {element_to_find}.")
    print(f"2. We search for this element in A.")
    print(f"3. The element {element_to_find} is found at index {i} in A.")
    print(f"\nTherefore, the rotation index is i = {i}.")

    # The problem asks to output each number in the final equation.
    # We will verify our finding and print the equation.
    print("\nVerification of the equation:")
    rotated_A = A[i:] + A[:i]
    print(f"A[{i}:] + A[:{i}] gives: {rotated_A}")
    print(f"B is:              {B}")
    # Using collections.deque for rotation is efficient but doesn't change overall complexity.
    # For instance, creating deque is O(n), rotation is O(k).
    # from collections import deque
    # A_deque = deque(A)
    # A_deque.rotate(-i)
    # assert list(A_deque) == B
    print(f"Is the equation satisfied? {rotated_A == B}")


# Example from the problem description
A = [1, 2, 3, 4, 5]
i_actual = 2
B = A[i_actual:] + A[:i_actual]

find_rotation_index(A, B)

# Another example
print("-" * 20)
A_2 = [30, 45, 62, 89, 10, 21]
i_actual_2 = 4
B_2 = A_2[i_actual_2:] + A_2[:i_actual_2]
find_rotation_index(A_2, B_2)