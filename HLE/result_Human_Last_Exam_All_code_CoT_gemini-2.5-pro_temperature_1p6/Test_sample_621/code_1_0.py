def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This function has a time complexity of O(n) due to the list.index() call,
    which performs a linear search. As argued in the analysis, this is optimal.
    
    Args:
        A: The original list of n unique integers.
        B: A list that is a cyclic shift of A.
    
    Returns:
        The integer rotation index i.
    """
    # Per the problem description, lists are non-empty and have the same length.
    # The first element of the rotated list B is our pivot.
    pivot_element = B[0]

    # The index of this pivot element in the original list A is the rotation index i.
    # The .index() method in Python performs a linear search, which is an O(n) operation.
    try:
        found_i = A.index(pivot_element)
    except ValueError:
        # This case should not occur if B is guaranteed to be a rotation of A.
        return -1

    return found_i

# --- Example Usage ---
# Let's define the lists A and B. B is A rotated by an index i.
A = [10, 20, 30, 40, 50, 60, 70]
i_actual = 4  # Let's use 4 as the true rotation index for our example
B = A[i_actual:] + A[:i_actual]

# Find the index using our algorithm
i_found = find_rotation_index(A, B)

# --- Output the results ---
print(f"Given List A: {A}")
print(f"Given List B: {B}")
print("-" * 30)

if i_found != -1:
    print(f"1. The first element of B is {B[0]}.")
    print(f"2. Searching for {B[0]} in A finds it at index {i_found}.")
    print(f"3. Therefore, the rotation index i = {i_found}.")
    print("\nVerification of the equation B = A[i:] + A[:i]:")
    print(f"  {B} = {A[i_found:]} + {A[:i_found]}")
else:
    print("Could not determine the rotation index.")
