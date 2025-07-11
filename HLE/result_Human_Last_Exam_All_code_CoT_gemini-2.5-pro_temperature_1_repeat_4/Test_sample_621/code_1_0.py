def find_rotation_index(A, B):
    """
    Finds the rotation index 'i' such that B is a cyclic shift of A.
    B = A[i:] + A[:i]

    This algorithm runs in O(n) time, which is optimal.

    Args:
        A (list): A list of n unique integers.
        B (list): A list that is a cyclic shift of A.
    """
    n = len(A)
    m = len(B)

    # Basic validation
    if n != m:
        print("Error: Lists A and B must have the same length.")
        return
    if n == 0:
        print("The lists are empty. The rotation index is 0.")
        return

    # Step 1: Find the candidate index 'i'.
    # Since all elements in A are unique, B[0] must match exactly one element in A.
    # We find the index 'i' of that element. This search takes O(n) time.
    try:
        first_element_of_B = B[0]
        candidate_i = A.index(first_element_of_B)
    except ValueError:
        print(f"Error: B is not a rotation of A, because {first_element_of_B} is not in A.")
        return

    # Step 2: Verify the rotation.
    # To be certain, we must check if the rest of the elements match the rotation.
    # This verification requires checking all n elements, taking O(n) time.
    # This step is why the problem has a lower bound of Omega(n).
    is_correct_rotation = True
    for j in range(n):
        if B[j] != A[(candidate_i + j) % n]:
            is_correct_rotation = False
            break

    # Step 3: Output the result.
    if is_correct_rotation:
        print(f"Given lists:")
        print(f"A = {A}")
        print(f"B = {B}")
        print("-" * 20)
        print(f"The rotation index is: i = {candidate_i}")
        
        # As requested, outputting the numbers in the final equation
        reconstructed_A_part1 = A[candidate_i:]
        reconstructed_A_part2 = A[:candidate_i]
        
        print(f"\nVerification of the equation B = A[i:] + A[:i]:")
        print(f"A[{candidate_i}:] = {reconstructed_A_part1}")
        print(f"A[:{candidate_i}] = {reconstructed_A_part2}")
        print(f"A[{candidate_i}:] + A[:{candidate_i}] = {reconstructed_A_part1 + reconstructed_A_part2}")
        print(f"This matches B, so the index is correct.")

    else:
        # This block should not be reached if B is guaranteed to be a rotation of A.
        print("Error: B is not a valid cyclic shift of A.")

# --- Example Usage ---
# Define the lists A and B. B is A rotated by i=3.
list_A = [10, 20, 30, 40, 50, 60]
list_B = [40, 50, 60, 10, 20, 30]

find_rotation_index(list_A, list_B)