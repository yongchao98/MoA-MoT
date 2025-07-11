def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    Args:
        A: The original list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i.
    """
    # The problem statement guarantees A and B are non-empty and B is a rotation of A.
    # If B is a left rotation of A by i, then B[0] must be equal to A[i].
    # B = [A[i], A[i+1], ..., A[n-1], A[0], ..., A[i-1]]
    # Since all elements in A are unique, finding the index of B[0] in A
    # gives us the rotation index i directly.
    
    # Get the first element of B. This takes O(1) time.
    target_element = B[0]

    # Find the index of this element in A. 
    # The list.index() method performs a linear search, which takes O(n) time
    # for an unsorted list.
    try:
        i = A.index(target_element)
    except ValueError:
        # This case is guaranteed not to happen by the problem statement.
        # If it did, it would mean B is not a rotation of A.
        return -1
        
    return i

# --- Example Usage ---

# Define the input lists A and B.
A = [10, 20, 30, 40, 50, 60]
# B is A shifted left by i=2.
# A[2:] = [30, 40, 50, 60]
# A[:2] = [10, 20]
# B = [30, 40, 50, 60, 10, 20]
B = [30, 40, 50, 60, 10, 20]

print(f"Given list A: {A}")
print(f"Given list B: {B}")
print("-" * 20)

# The core of the algorithm is to solve for 'i' in the equation B[0] = A[i]
i = find_rotation_index(A, B)

# Per instructions to output numbers in the final equation:
print("The algorithm is based on the equation: B[0] = A[i]")
print(f"In this example, B[0] = {B[0]}.")
print(f"We search for the element {B[0]} in list A.")
print(f"The index of {B[0]} in A is {i}.")
print(f"Therefore, the rotation index is i = {i}.")
print("-" * 20)
print(f"Final Equation: {B[0]} = A[{i}]")

<<<A>>>