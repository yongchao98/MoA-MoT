def find_rotation_index(A, B):
    """
    Finds the rotation index i for two lists A and B, where B is a rotation of A.
    
    The problem asks for the best possible time complexity. The relationship
    B = A[i:] + A[:i] implies that the first element of B, B[0], must be A[i].
    Since all integers in A are unique, finding the index of B[0] in A will
    uniquely determine i.
    
    This task is equivalent to searching for an element in an unsorted list,
    which has a worst-case time complexity lower bound of Omega(n). An O(n)
    algorithm is therefore optimal. This implementation achieves O(n) by using
    the list.index() method.
    
    Args:
        A (list): The original list of n unique integers.
        B (list): A list created by rotating A, i.e., B = A[i:] + A[:i].
    """
    if not A:
        print("Input list A cannot be empty.")
        return

    # The rotation index 'i' is the index of the first element of B in list A.
    # This search operation is O(n).
    start_element = B[0]
    try:
        i = A.index(start_element)
    except ValueError:
        # This case won't happen based on the problem description but is good practice.
        print(f"Error: B is not a valid rotation of A, as {start_element} is not in A.")
        return

    print(f"Given A = {A}")
    print(f"Given B = {B}")
    print(f"The rotation index is i = {i}.")
    print("\nThis means B = A[i:] + A[:i]. Let's verify:")
    
    # To fulfill the request to "output each number in the final equation",
    # we construct string representations of the lists.
    b_str = ", ".join(map(str, B))
    a_part1_str = ", ".join(map(str, A[i:]))
    a_part2_str = ", ".join(map(str, A[:i]))
    
    print(f"[{b_str}] = [{a_part1_str}] + [{a_part2_str}]")


# --- Example Usage ---
# Define the lists A and B for demonstration.
A = [10, 20, 30, 40, 50, 60, 70]
# B is A rotated by i=3, so B = A[3:] + A[:3]
B = [40, 50, 60, 70, 10, 20, 30]

find_rotation_index(A, B)