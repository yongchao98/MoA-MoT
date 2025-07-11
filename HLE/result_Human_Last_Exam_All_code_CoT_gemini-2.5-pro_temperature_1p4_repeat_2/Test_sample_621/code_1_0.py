def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    The algorithm has a time complexity of O(n) due to the list search.

    Args:
        A: A list of n unique integers.
        B: A list that is a rotation of A.
        
    Returns:
        The integer rotation index i, or -1 if inputs are invalid.
    """
    # If lists are empty or have different lengths, it's not a valid rotation.
    if not A or not B or len(A) != len(B):
        # The problem statement guarantees valid inputs, but robust code should check.
        return -1
        
    # The first element of B must be A[i].
    first_element_of_B = B[0]
    
    # Since all elements in A are unique, we just need to find the index
    # of B[0] in A. list.index() performs a linear search, which is O(n).
    try:
        i = A.index(first_element_of_B)
    except ValueError:
        # This case won't happen if B is guaranteed to be a rotation of A.
        return -1
        
    return i

# --- Example Usage ---
# Define the input lists A and B
A = [3, 8, 1, 9, 4, 12]
B = [9, 4, 12, 3, 8, 1]

# Find the rotation index using the function
i = find_rotation_index(A, B)

# --- Outputting the results as requested ---
print(f"Given list A: {A}")
print(f"Given rotated list B: {B}")

if i != -1:
    print(f"\nThe rotation index is i = {i}.")
    print("The equation is B = A[i:] + A[:i]")
    
    # To show the equation with the numbers
    part1 = A[i:]
    part2 = A[:i]
    print(f"Substituting the values: {B} = {part1} + {part2}")
else:
    print("Could not find a valid rotation index.")
