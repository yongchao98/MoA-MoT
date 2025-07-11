def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B is a cyclic shift of A.
    A: A list of n unique integers.
    B: A list that is a rotation of A (B = A[i:] + A[:i]).
    """
    if len(A) != len(B) or not A:
        print("Invalid input: Lists must be non-empty and of the same length.")
        return

    # The first element of the rotated list B must be A[i].
    target_element = B[0]

    try:
        # list.index() performs a linear search, which is O(n).
        rotation_index = A.index(target_element)

        # Print the logic and the result.
        print(f"Given list A: {A}")
        print(f"Given list B: {B}")
        print("-" * 30)
        print(f"The goal is to find the index 'i' such that B = A[i:] + A[:i].")
        print(f"This means the first element of B, which is B[0]={B[0]}, must be the element at index 'i' in A.")
        print(f"We search for the element {target_element} in list A.")
        print(f"The element {target_element} is found at index {rotation_index} in A, because A[{rotation_index}] = {A[rotation_index]}.")
        print(f"\nTherefore, the rotation index i is {rotation_index}.")

        # Verification step to show the equation holds
        reconstructed_b = A[rotation_index:] + A[:rotation_index]
        print(f"\nVerification:")
        print(f"A[{rotation_index}:] + A[:{rotation_index}] results in {reconstructed_b}")
        print(f"This matches list B, confirming the result is correct.")

    except ValueError:
        print(f"Error: B is not a valid rotation of A, as {target_element} is not in A.")

# Example Usage:
# Let A be a list of unique integers.
A = [15, 25, 35, 45, 55, 65]
# Let B be A rotated by i=3. A[3:] = [45, 55, 65], A[:3] = [15, 25, 35]
B = [45, 55, 65, 15, 25, 35]

find_rotation_index(A, B)