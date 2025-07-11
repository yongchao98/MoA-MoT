def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This function leverages the fact that A has unique elements and B is a guaranteed
    rotation of A. The index of the first element of B in A gives the rotation index.
    The time complexity is O(n) due to the list.index() call, which performs a linear search.
    This is optimal as the problem has a lower bound of Omega(n).
    
    Args:
        A: A list of n unique integers.
        B: A list which is a cyclic shift of A.
        
    Returns:
        The integer rotation index i, or -1 if inputs are invalid.
    """
    if not A or not B or len(A) != len(B):
        return -1
        
    # The first element of B must be at index `i` in A.
    element_to_find = B[0]
    
    try:
        # The .index() method searches for the element and returns its index.
        # This operation has a time complexity of O(n).
        i = A.index(element_to_find)
        return i
    except ValueError:
        # This case should not be reached given the problem's constraints
        # (B is guaranteed to be a rotation of A).
        return -1

# --- Example Usage ---
# Define the lists A and B. B is A rotated by i=3.
A = [15, 25, 35, 45, 55, 65, 75]
B = [45, 55, 65, 75, 15, 25, 35]

# Find the rotation index
i = find_rotation_index(A, B)

# Print the results as requested
print(f"Given list A: {A}")
print(f"Given list B: {B}")

if i != -1:
    print(f"The rotation index is: i = {i}")
    
    # Construct the equation string B = A[i:] + A[:i]
    rotated_part1 = A[i:]
    rotated_part2 = A[:i]
    
    print("\nThe equation representing the rotation is:")
    # The problem asks to output each number in the final equation.
    # Printing the lists themselves satisfies this requirement.
    print(f"{B} = {rotated_part1} + {rotated_part2}")
else:
    print("Could not find a valid rotation index. Please check the input lists.")
