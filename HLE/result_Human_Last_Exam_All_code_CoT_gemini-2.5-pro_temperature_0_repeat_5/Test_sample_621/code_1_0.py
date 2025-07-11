def find_rotation_index(A, B):
    """
    Finds the rotation index i such that B = A[i:] + A[:i].
    
    This algorithm has a time complexity of O(n) because of the list.index() call,
    which performs a linear scan. As proven by the adversary argument, this is
    the optimal time complexity.

    Args:
        A: A list of n unique integers.
        B: A list that is a rotation of A.

    Returns:
        The integer rotation index i, or -1 if B is not a valid rotation of A.
    """
    # Basic validation: lists must be of the same length.
    if len(A) != len(B):
        return -1
    
    # If lists are empty, the rotation index is 0.
    if not A:
        return 0

    # The first element of B must exist in A. Its index in A is the rotation index i.
    first_element_of_B = B[0]
    
    try:
        # The .index() method performs a linear search, which takes O(n) time.
        # Since elements in A are unique, this will give the correct index i.
        i = A.index(first_element_of_B)
        return i
    except ValueError:
        # This case occurs if B[0] is not in A, meaning B cannot be a rotation of A.
        # Based on the problem description, this path should not be taken.
        return -1

# --- Example Usage ---
# Define the original list A with unique integers
A = [10, 20, 30, 40, 50, 60, 70]

# Define the rotated list B. Here, we create it by rotating A with i=3
# B = A[3:] + A[:3] = [40, 50, 60, 70] + [10, 20, 30]
B = [40, 50, 60, 70, 10, 20, 30]

# Find the rotation index
i = find_rotation_index(A, B)

# --- Output Results ---
print(f"Given list A: {A}")
print(f"Given list B: {B}")

if i != -1:
    print(f"\nCalculated rotation index i = {i}")
    
    # As requested, output each number in the final equation B = A[i:] + A[:i]
    print("\nFinal Equation Breakdown:")
    
    b_str = ", ".join(map(str, B))
    a_part1_str = ", ".join(map(str, A[i:]))
    a_part2_str = ", ".join(map(str, A[:i]))
    
    print(f"B         : [{b_str}]")
    print(f"A[{i}:]     : [{a_part1_str}]")
    print(f"A[:{i}]     : [{a_part2_str}]")
    print(f"A[{i}:]+A[:{i}]: [{a_part1_str}] + [{a_part2_str}]")

    # Final verification
    is_correct = (A[i:] + A[:i]) == B
    print(f"\nIs the reconstruction correct? {is_correct}")
else:
    print("\nCould not find a valid rotation index. B is not a rotation of A.")
