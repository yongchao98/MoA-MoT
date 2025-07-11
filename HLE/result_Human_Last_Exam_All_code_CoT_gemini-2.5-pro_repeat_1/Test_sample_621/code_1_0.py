def find_rotation_index(A, B):
    """
    Given a list A and a list B, where B is a cyclic shift of A (B = A[i:] + A[:i]),
    this function finds the rotation index i.
    The lists are assumed to contain unique integers.

    Args:
      A: The original list of unique integers.
      B: The rotated list.

    Returns:
      The integer rotation index i.
    """
    n = len(A)
    # Handle edge cases
    if n == 0:
        return 0
    if n != len(B):
        raise ValueError("Input lists A and B must have the same length.")

    # To find the index i efficiently, we can create a mapping from
    # each element in A to its index. A hash map (dictionary) is perfect for this.
    # This step has a time complexity of O(n).
    a_element_to_index_map = {element: index for index, element in enumerate(A)}

    # The first element of the rotated list B, B[0], must be the element A[i].
    # So, we can find i by looking up the index of B[0] in our map of A.
    # This lookup is O(1) on average.
    first_element_in_B = B[0]

    # Find the index in the original list A
    rotation_index = a_element_to_index_map.get(first_element_in_B)

    if rotation_index is None:
        # This case should not be reached if B is guaranteed to be a rotation of A.
        raise ValueError("B is not a valid rotation of A.")

    return rotation_index

# --- Example Usage ---
# Let A be a list of numbers.
A = [10, 20, 30, 40, 50, 60]
# Let B be A rotated by i=3.
# A[3:] = [40, 50, 60]
# A[:3] = [10, 20, 30]
# B = [40, 50, 60, 10, 20, 30]
B = [40, 50, 60, 10, 20, 30]

# We call our function to find the rotation index.
i = find_rotation_index(A, B)

# The prompt requires printing the numbers in the final equation.
# The "equation" is B = A[i:] + A[:i], and the number we solved for is i.
# So, we will print the value of i.
print(i)