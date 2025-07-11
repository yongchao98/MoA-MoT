def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].

    The time complexity of this function is O(n) because of the `A.index(target_element)`
    call, which performs a linear search on list A. As explained in the reasoning,
    this is the best possible time complexity because, in the worst case, any
    algorithm must inspect all elements of A to locate the position of B[0].

    Args:
        A (list): The original list of n unique integers.
        B (list): A list created by cyclically shifting A.

    Returns:
        int: The rotation index i, or -1 if inputs are invalid.
    """
    n = len(A)
    if n != len(B) or n == 0:
        return -1 # Invalid input

    # The first element of the rotated list B must be A[i].
    target_element = B[0]

    try:
        # Find the index of the target element in A.
        # This is a linear search with O(n) time complexity.
        i = A.index(target_element)
        return i
    except ValueError:
        # This case should not happen if B is guaranteed to be a rotation of A.
        return -1

# --- Example Usage ---

# Define the original list A
A = [10, 20, 30, 40, 50, 60, 70, 80]
# Define the rotation index we want to find
original_i = 5
# Create the rotated list B = A[5:] + A[:5]
B = A[original_i:] + A[:original_i]

print(f"Original list A: {A}")
print(f"Rotated list B:  {B}")
print("-" * 20)

# Use the function to find the rotation index
found_i = find_rotation_index(A, B)

if found_i != -1:
    print(f"The algorithm found the rotation index i = {found_i}")

    # Verification as per the prompt's "final equation" hint
    # B = A[i:] + A[:i]
    print("\nVerifying the result...")
    part1 = A[found_i:]
    part2 = A[:found_i]
    reconstructed_B = part1 + part2

    print(f"The first part of the equation is A[i:] = A[{found_i}:] = {part1}")
    print(f"The second part of the equation is A[:i] = A[:{found_i}] = {part2}")
    print(f"The reconstructed list B = {part1} + {part2} = {reconstructed_B}")
    print(f"This matches the input B: {B == reconstructed_B}")
else:
    print("Could not find the rotation index or inputs were invalid.")
