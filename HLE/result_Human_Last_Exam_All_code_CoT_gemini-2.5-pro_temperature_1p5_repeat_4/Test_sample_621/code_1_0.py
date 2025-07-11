def find_rotation_index(A, B):
    """
    Finds the rotation index i for two lists A and B, where B is a cyclic shift of A.
    It is guaranteed that A has n unique integers and B is a rotation of A.

    Args:
        A (list): The original list of unique integers.
        B (list): The rotated list.

    Returns:
        int: The rotation index i, or -1 if inputs are invalid.
    """
    if not A or len(A) != len(B):
        print("Invalid input: Lists must be non-empty and of the same length.")
        return -1
        
    n = len(A)
    
    # 1. The first element of B, B[0], corresponds to A[i].
    value_to_find = B[0]
    
    # 2. Find the index of this element in A. Since all elements in A are unique,
    #    this index will also be unique. Python's list.index() does this search.
    #    This search operation is O(n) as it may need to scan the entire list.
    try:
        i = A.index(value_to_find)
    except ValueError:
        print(f"Error: The first element of B ({value_to_find}) was not found in A.")
        return -1

    # 3. Print the results and the final equation.
    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print("\nStep 1: Identify the element in A that starts the rotated sequence B.")
    print(f"This element is B[0], which is {value_to_find}.")
    print("\nStep 2: Find the index of this element in A.")
    print(f"The value {value_to_find} is found at index {i} in list A.")
    print("\nThis means the rotation index is i =", i)
    
    print("\nThe final equation is B = A[i:] + A[:i]")
    # We will print the numbers in the equation for clarity
    a_part1_str = ", ".join(map(str, A[i:]))
    a_part2_str = ", ".join(map(str, A[:i]))
    b_str = ", ".join(map(str, B))

    print(f"[{b_str}] = [{a_part1_str}] + [{a_part2_str}]")

    # The logic confirms that the complexity is O(n), and it cannot be improved.
    
    return i

# --- Example Usage ---
# Let's define two lists A and B where B is A rotated by i=3.
A = [10, 20, 30, 40, 50, 60]
i_actual = 3
B = A[i_actual:] + A[:i_actual]

# Run the algorithm to find i
find_rotation_index(A, B)
