def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    This algorithm is O(n) because list.index() performs a linear search.
    This is the optimal time complexity, as finding an element in an unsorted
    list has a lower bound of Omega(n).

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The integer rotation index i.
    """
    if not A or not B or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same size.")
        return -1

    # Take the first element of B as the pivot.
    pivot = B[0]

    try:
        # Find the index of the pivot element in A.
        # This is the O(n) step.
        i = A.index(pivot)
        return i
    except ValueError:
        # This case should not be reached if B is truly a rotation of A.
        print(f"Error: Pivot element {pivot} from B not found in A.")
        return -1

# --- Example Usage ---
A = [10, 20, 30, 40, 50, 60]
i_true = 4
B = A[i_true:] + A[:i_true]

print(f"Given list A: {A}")
print(f"Given list B: {B}")
print("-" * 20)

# Find the rotation index
found_i = find_rotation_index(A, B)

if found_i != -1:
    print(f"Found rotation index i = {found_i}")

    # The prompt requests to output the numbers in the final equation.
    # Let's verify and print the equation B = A[i:] + A[:i].
    print("\nVerifying the equation:")
    # Print each part of the equation with its full list content.
    print(f"B           = {B}")
    rotated_A = A[found_i:] + A[:found_i]
    print(f"A[{found_i}:]+A[:{found_i}] = {A[found_i:]} + {A[:found_i]}")
    print(f"Result      = {rotated_A}")
    
    assert B == rotated_A
    print("\nVerification successful.")
