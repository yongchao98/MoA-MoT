def find_rotation_index(A, B):
    """
    Given a list A of n unique integers and a rotated list B = A[i:] + A[:i],
    this function finds the rotation index i.

    The algorithm has a time complexity of O(n) because it involves a linear
    search in the unsorted list A. As explained in the reasoning, this is
    the optimal time complexity for this problem.
    """
    n = len(A)
    if n == 0:
        print("Lists are empty, cannot determine rotation index.")
        return

    if n != len(B):
        print("Lists have different lengths, B cannot be a rotation of A.")
        return

    # The first element of B must be the element at index 'i' in A.
    target_element = B[0]

    # We find the index of this target element in A.
    # The list.index() method performs a linear search, which is O(n).
    try:
        i = A.index(target_element)

        # Output the results as requested
        print(f"Input list A: {A}")
        print(f"Input list B: {B}")
        print("-" * 20)
        print(f"To find the rotation index 'i', we take the first element of B, which is {target_element}.")
        print(f"We search for {target_element} in list A.")
        print(f"The element {target_element} is found at index {i} in A.")
        print(f"\nTherefore, the rotation index is i = {i}.")
        
        # Verification and printing the equation with numbers
        print("\nVerification of the equation B = A[i:] + A[:i]:")
        a_suffix = A[i:]
        a_prefix = A[:i]
        
        # Printing each number in the final equation
        print(f"B = {B}")
        print(f"A[i:] = A[{i}:] = {a_suffix}")
        print(f"A[:i] = A[:{i}] = {a_prefix}")
        print(f"A[i:] + A[:i] = {a_suffix} + {a_prefix} = {a_suffix + a_prefix}")
        
        if (a_suffix + a_prefix) == B:
            print("\nThe equation holds true.")
        else:
            # This part should not be reached if B is guaranteed to be a rotation of A
            print("\nError: The reconstructed list does not match B.")

    except ValueError:
        # This error occurs if B[0] is not in A, meaning B is not a rotation of A.
        print(f"Error: The element {target_element} from B was not found in A.")
        print("B is not a valid rotation of A.")

# --- Example Usage ---
# Define the original list A
A = [15, 25, 35, 45, 55, 65, 75]

# Define the rotation index i
i = 4

# Create the rotated list B based on A and i
B = A[i:] + A[:i]

# Call the function to find i
find_rotation_index(A, B)