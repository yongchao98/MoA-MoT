def find_rotation_index(A, B):
    """
    Given a list A of n unique integers and a list B which is a cyclic shift of A,
    this function finds the index i such that B = A[i:] + A[:i].

    The algorithm runs in O(n) time.
    """
    n = len(A)

    # Basic validation
    if n != len(B):
        print("Error: Lists A and B must have the same length.")
        return
    
    if n == 0:
        # Handle the edge case of empty lists
        print("For empty lists, the shift index is 0.")
        print("Equation: [] = [] + []")
        return
    
    # Since all elements in A are unique, the first element of B must match
    # exactly one element in A. Let's find its index.
    # This is our candidate for the shift index 'i'.
    # This step has a time complexity of O(n).
    first_element_of_B = B[0]
    try:
        i = A.index(first_element_of_B)
    except ValueError:
        print(f"Error: The first element of B ({first_element_of_B}) was not found in A.")
        print("Therefore, B cannot be a cyclic shift of A.")
        return

    # To be certain, we must verify that the rest of the list matches the shift.
    # The problem statement guarantees that B is a shift of A, so this check
    # is technically for correctness and demonstration.
    # This step of slicing, concatenation, and comparison also takes O(n).
    expected_B = A[i:] + A[:i]
    if B == expected_B:
        print(f"The list A is: {A}")
        print(f"The list B is: {B}")
        print(f"The rotation index 'i' is: {i}")
        print("\nFinal Equation:")
        # The prompt requires printing each number in the final equation.
        # The f-string representation of a list accomplishes this.
        print(f"{B} = {A[i:]} + {A[:i]}")
        print(f"Which corresponds to the formula: B = A[{i}:] + A[:{i}]")
    else:
        # This block should not be reached if the input is valid.
        print(f"Verification failed. B is not the expected rotation of A for index i={i}.")

# --- Example Usage ---
# Let's define two lists A and B where B is a cyclic shift of A.
A = [10, 20, 30, 40, 50, 60]
# We will create B by shifting A by i=2
# B = A[2:] + A[:2] = [30, 40, 50, 60] + [10, 20]
B = [30, 40, 50, 60, 10, 20]

# Run the algorithm to find the index 'i'
find_rotation_index(A, B)