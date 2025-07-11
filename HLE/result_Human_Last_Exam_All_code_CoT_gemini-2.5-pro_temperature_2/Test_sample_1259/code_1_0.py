import math
import sys

def solve():
    """
    Solves the problem of restoring sequence A from its set of subarray GCDs S.
    """
    try:
        # Read the number of elements in S (not strictly needed if reading the whole line)
        if sys.stdin.isatty():
            # Handle potential manual input format difference
            _ = int(sys.stdin.readline())
        
        # Read the elements of S from a single line
        s_str = sys.stdin.readline().strip().split()
        if not s_str:
            print("-1")
            return

        s_set = set(map(int, s_str))
        s_list = sorted(list(s_set))
        m = len(s_list)

    except (IOError, ValueError):
        print("-1")
        return

    # The minimum element of S must be the gcd of the whole sequence A.
    g = s_list[0]

    # This minimum element 'g' must divide all elements of S.
    # This is a necessary and sufficient condition.
    # Check if this condition holds.
    for x in s_list:
        if x % g != 0:
            print("-1")
            return

    # If the condition holds, a solution can be constructed.
    # A simple construction is to interleave g with the other elements of S.
    # If S = {g}, then A = [g].
    if m == 1:
        print(1)
        print(g)
    else:
        # For |S| > 1, we construct A = [s'_1, g, s'_2, g, ...].
        # where s'_i are elements of S other than g.
        result_a = []
        for i in range(1, m):
            result_a.append(s_list[i])
            result_a.append(g)
        
        # Output the length of the constructed sequence A.
        print(len(result_a))
        # Output the elements of A, separated by spaces.
        print(*result_a)

solve()
<<<H>>>