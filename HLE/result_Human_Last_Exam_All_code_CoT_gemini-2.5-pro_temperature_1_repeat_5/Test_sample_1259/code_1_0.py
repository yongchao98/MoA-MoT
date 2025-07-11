import math

def solve_gcd_construction(S):
    """
    Given a list of integers S, this function attempts to restore a sequence A
    such that the set of GCDs of all subarrays of A is equal to S.

    The logic is based on the finding that a solution exists if and only if
    the minimum element of S divides all other elements of S.

    If a solution exists, a valid sequence A can be constructed by interleaving
    the minimum element `g` with the other elements of S. For example, if
    S = {g, s2, s3}, a valid A would be [g, s2, g, s3].

    The function prints the constructed sequence A as specified. The instruction
    "output each number in the final equation" is interpreted as printing the
    elements of the sequence A.
    """
    if not S:
        print("0")
        print()
        return

    # Use a set to remove duplicates and then sort to easily find the minimum.
    s_list = sorted(list(set(S)))
    m = len(s_list)

    if m == 1:
        # If there's only one element, A is just that element.
        print("1")
        print(s_list[0])
        return

    g = s_list[0]

    # Check if the necessary condition holds (min element divides all others).
    # If not, no solution exists.
    for x in s_list:
        if x % g != 0:
            print("-1")
            return

    # Construct the sequence A by interleaving g with other elements.
    # A = [g, s_2, g, s_3, ..., g, s_k]
    result_a = []
    for i in range(m):
        result_a.append(s_list[i])
        if i > 0: # Append g after every element except the first one (which is g itself)
             result_a.append(g)

    # My explanation describes A=[g, s2, g, s3, ...]. Let's implement that.
    # It's simpler and also correct.
    result_a_interleaved = []
    result_a_interleaved.append(g)
    for i in range(1, m):
        result_a_interleaved.append(s_list[i])
        result_a_interleaved.append(g)


    print(len(result_a_interleaved))
    print(*result_a_interleaved)


# Example usage:
# S = [2, 4, 6, 12]
# solve_gcd_construction(S)
# Expected Output:
# 7
# 2 4 2 6 2 12 2

# S = [1, 6, 10, 15]
# solve_gcd_construction(S)
# Expected Output:
# 7
# 1 6 1 10 1 15 1

# S = {6, 12, 18, 36}
# solve_gcd_construction(S)
# Expected Output:
# 7
# 6 12 6 18 6 36 6

# Note: The code is provided in a function for clarity. 
# For a typical competitive programming setup, you would read from stdin.
# For example:
# m_input = int(input())
# s_input = list(map(int, input().split()))
# solve_gcd_construction(s_input)