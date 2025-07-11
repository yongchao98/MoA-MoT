import random

def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    This is an O(n) implementation, which is optimal.

    Args:
        A: A list of n unique integers.
        B: A list which is a cyclic shift of A.

    Returns:
        The integer rotation index i.
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: lists must be non-empty and have the same size.")
        return None

    # The first element of B must be A[i].
    target_value = B[0]

    try:
        # list.index() performs a linear search, which is O(n).
        i = A.index(target_value)
    except ValueError:
        print(f"Error: B is not a rotation of A, as {target_value} is not in A.")
        return None

    # As per the problem description, we can assume B is always a rotation of A.
    # The 'try...except' block is for robustness in a general case.
    
    print(f"Given lists:")
    print(f"A = {A}")
    print(f"B = {B}")
    print("-" * 20)
    
    # The problem is solved by finding i such that A[i] = B[0].
    print("The final equation we need to solve is: A[i] = B[0]")
    
    print("Each number in the final equation:")
    print(f"B[0] = {B[0]}")
    print(f"The index 'i' where this value is found in A is: {i}")
    print(f"The value at A[i] (which is A[{i}]) is: {A[i]}")
    print("-" * 20)
    
    print(f"The rotation index is i = {i}")
    return i

# --- Example Usage ---
# Let's create a sample A and a rotated version B
n = 10
# Generate a list of n unique integers
A_list = random.sample(range(1, 101), n) 

# Generate a random rotation index i
rotation_i = random.randint(0, n - 1)

# Create the rotated list B
B_list = A_list[rotation_i:] + A_list[:rotation_i]

# Find the rotation index using the function
find_rotation_index(A_list, B_list)