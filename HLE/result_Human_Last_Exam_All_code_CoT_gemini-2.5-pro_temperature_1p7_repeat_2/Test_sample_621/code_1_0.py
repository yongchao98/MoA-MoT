def find_rotation_index(A, B):
    """
    Given a list A of unique integers and a list B which is a cyclic shift of A,
    this function finds the rotation index i such that B = A[i:] + A[:i].

    The algorithm has a time complexity of O(n) because of the list.index() search,
    which is proven to be optimal.
    """
    n = len(A)
    m = len(B)

    # Basic validation
    if n != m:
        print("Error: Lists must be of the same length.")
        return

    if n == 0:
        # Handle the edge case of empty lists.
        print("Given A = [] and B = []")
        print("The rotation index 'i' is 0.")
        print("[] = [][0:] + [][:0]")
        return

    try:
        # The first element of B must correspond to A[i].
        # B[0] == A[i]
        first_element_of_B = B[0]
        
        # Find the index of B[0] in A. This takes O(n) time.
        # Because all elements in A are unique, this gives the correct i.
        i = A.index(first_element_of_B)
        
        # As per the problem statement, B is guaranteed to be a rotation of A,
        # so no further verification is strictly needed.

        print(f"Given A = {A} and B = {B}")
        print(f"The rotation index 'i' is {i}.")
        print("The final equation is B = A[i:] + A[:i], which is:")
        # Output the equation with the actual numbers and slices
        print(f"{B} = {A[i:]} + {A[:i]}")

    except ValueError:
        # This case happens if the first element of B is not in A.
        print(f"Error: The element {B[0]} from B was not found in A.")
        print("B cannot be a cyclic shift of A.")


# --- Example Usage ---
# Let's define two lists where B is A rotated by i=3
list_A = [15, 25, 35, 45, 55, 65]
# When i=3, A[3:] = [45, 55, 65] and A[:3] = [15, 25, 35]
# B = [45, 55, 65, 15, 25, 35]
list_B = [45, 55, 65, 15, 25, 35]

find_rotation_index(list_A, list_B)