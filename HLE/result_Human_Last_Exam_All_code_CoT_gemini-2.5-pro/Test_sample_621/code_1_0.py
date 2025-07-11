def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A: A list of n unique integers.
    B: A list that is a cyclic shift of A.
    Returns the rotation index i.
    """
    if not A or not B or len(A) != len(B):
        return -1  # Or raise an error for invalid input

    # The first element of the rotated list B must be A[i].
    first_element_of_B = B[0]

    # Find the index of this element in A.
    # The .index() method performs a linear search, which is O(n).
    try:
        rotation_index = A.index(first_element_of_B)
        return rotation_index
    except ValueError:
        # This case should not happen based on the problem description,
        # as B is guaranteed to be a rotation of A.
        return -1

# --- Example Usage ---
# Define the original list A and the rotated list B
A = [10, 20, 30, 40, 50, 60]
# Let's choose a rotation index, for example i = 2
i_actual = 2
B = A[i_actual:] + A[:i_actual]

# Find the rotation index using the algorithm
i_found = find_rotation_index(A, B)

# Print the results and the final equation
print(f"Original list A: {A}")
print(f"Rotated list B:  {B}")
print(f"Found rotation index i = {i_found}")

# Verify and print the equation with the numbers
if i_found != -1:
    reconstructed_B = A[i_found:] + A[:i_found]
    # The problem asks to output each number in the final equation
    print("\nVerification using the found index:")
    print(f"{B} == {A[i_found:]} + {A[:i_found]}")
    # Check if the reconstruction is correct
    if B == reconstructed_B:
        print("The equation holds true.")
    else:
        print("Error: The reconstruction does not match B.")