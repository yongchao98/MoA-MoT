import random

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i, or -1 if B is not a rotation of A
        or if lists are empty.
    """
    if not A or not B or len(A) != len(B):
        return -1

    # The first element of B must be some element A[i].
    # Since all elements in A are unique, we can find this i
    # by searching for B[0] in A.
    first_element_of_B = B[0]

    try:
        # list.index() performs a linear search, which is O(n).
        # This is the dominant operation.
        i = A.index(first_element_of_B)
        
        # As per the problem description, B is guaranteed to be a rotation of A.
        # If we wanted to be defensive, we could verify the rotation:
        # assert B == A[i:] + A[:i]
        
        return i
    except ValueError:
        # This case should not happen if B is guaranteed to be a rotation of A.
        return -1

# --- Example Usage ---

# Define an initial list A
A = [10, 20, 30, 40, 50, 60]
n = len(A)

# Pick a random rotation index i
# For demonstration, let's use a fixed index
# i = random.randint(0, n - 1)
i_original = 4 

# Create the rotated list B
B = A[i_original:] + A[:i_original]

print(f"Given list A: {A}")
print(f"Given list B: {B}")
print("-" * 20)

# Find the rotation index using the algorithm
found_i = find_rotation_index(A, B)

if found_i != -1:
    print(f"The algorithm found the rotation index: i = {found_i}")
    
    # Outputting the numbers in the final equation as verification
    print("\nVerification of the rotation equation B = A[i:] + A[:i]:")
    print(f"B = {B}")
    reconstructed_B = A[found_i:] + A[:found_i]
    # This line shows the equation with all numbers from the lists
    print(f"A[{found_i}:] + A[:{found_i}] = {A[found_i:]} + {A[:found_i]} = {reconstructed_B}")

    if B == reconstructed_B:
        print("\nVerification successful: The lists match.")
    else:
        print("\nVerification failed.")
else:
    print("Could not find the rotation index. B may not be a rotation of A.")
