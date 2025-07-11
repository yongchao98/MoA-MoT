def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    The algorithm has a time complexity of O(n) due to the list search.

    Args:
        A: A list of n unique integers.
        B: A list that is a cyclic shift of A.

    Returns:
        The rotation index i.
    """
    if not A or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same length.")
        return -1
    
    n = len(A)
    if n == 0:
        return 0

    # The core idea: The first element of the rotated list B must be the
    # element at index 'i' in the original list A.
    first_element_of_B = B[0]
    
    print(f"List A: {A}")
    print(f"List B: {B}")
    print("-" * 20)
    print(f"Step 1: Get the first element of B, which is {first_element_of_B}.")

    try:
        # Step 2: Find the index of this element in A.
        # This list.index() call is the O(n) operation.
        i = A.index(first_element_of_B)
        
        print(f"Step 2: Find the index of {first_element_of_B} in A. It is at index {i}.")
        print(f"Conclusion: The rotation index is {i}.")
        
        # The problem asks to output the numbers in the final equation.
        # We will show the verification B == A[i:] + A[:i] with the actual list values.
        print("\n--- Final Equation Verification ---")
        rotated_part1 = A[i:]
        rotated_part2 = A[:i]
        reconstructed_B = rotated_part1 + rotated_part2
        
        print(f"Original B:    {B}")
        print(f"A[{i}:] + A[:{i}]: {rotated_part1} + {rotated_part2}")
        print(f"Reconstructed B: {reconstructed_B}")

        # This check should always pass given the problem statement.
        if reconstructed_B == B:
            print("\nVerification successful.")
        else:
            print("\nVerification failed. This should not happen.")
            
        return i
        
    except ValueError:
        # This error occurs if B[0] is not in A, which contradicts the problem statement.
        print(f"Error: The element {first_element_of_B} from B was not found in A.")
        return -1

# --- Example ---
# Define an original list A and a rotation index i
A_list = [10, 20, 30, 40, 50, 60]
rotation_i = 3

# Create the rotated list B
B_list = A_list[rotation_i:] + A_list[:rotation_i]

# Find the rotation index using the function
find_rotation_index(A_list, B_list)