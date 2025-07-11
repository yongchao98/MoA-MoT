import random

def find_rotation_index(A, B):
    """
    Given a list A of n unique integers and a list B which is a cyclic shift of A,
    this function finds the rotation index i such that B = A[i:] + A[:i].

    The time complexity of this function is O(n) because of the list.index() method,
    which performs a linear search. As explained in the analysis, this is the
    optimal time complexity for this problem.

    Args:
        A: A list of n unique integers.
        B: A list which is a cyclic shift of A.
    """
    print(f"Original list A: {A}")
    print(f"Rotated list B: {B}")
    print("-" * 20)

    if not A:
        print("Lists are empty, rotation index is 0.")
        return

    # From the definition of a cyclic shift, we know that B[0] must be equal to A[i].
    # So, we need to find the index 'i' in A where the element is equal to B[0].
    target_value = B[0]

    print(f"The key relationship is B[0] = A[i].")
    print(f"Here, B[0] = {target_value}.")
    print(f"We need to find the index 'i' in A such that A[i] = {target_value}.")

    # The list.index() method performs a linear scan of the list to find the
    # index of the first occurrence of an element. This takes O(n) time.
    try:
        i = A.index(target_value)
        print(f"Searching for {target_value} in A...")
        print(f"Found {target_value} at index {i} in A, since A[{i}] = {A[i]}.")
        print(f"Thus, the final equation is {B[0]} = A[{i}].")
        print(f"The rotation index is i = {i}.")
    except ValueError:
        # This part of the code should not be reached if the input is valid as per the problem description.
        print(f"Error: The element {target_value} (from B[0]) was not found in A.")


# --- Example Usage ---
# You can modify these lists to test with your own examples.
# Create a sample list A with n unique integers
n = 10
A = random.sample(range(1, 101), n)

# Choose a random rotation index i
i = random.randint(0, n - 1)

# Create B by rotating A by i
B = A[i:] + A[:i]

find_rotation_index(A, B)