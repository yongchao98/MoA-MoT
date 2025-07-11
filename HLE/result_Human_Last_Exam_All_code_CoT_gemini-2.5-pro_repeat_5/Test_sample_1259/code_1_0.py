import math
from functools import reduce

def solve():
    """
    Solves the problem of restoring sequence A from its GCD set S.
    The solution is based on the finding that a valid sequence A exists if and only if
    the minimum element of S is equal to the GCD of all elements in S.
    """
    # It is assumed the input is a single line of space-separated integers representing the set S.
    try:
        s_list_input = list(map(int, input().split()))
    except (IOError, ValueError):
        # Handle cases with no input or invalid input format.
        print("-1")
        return

    if not s_list_input:
        print("-1")
        return

    # A set has unique elements. We sort the list to easily find the minimum.
    s_list = sorted(list(set(s_list_input)))
    m = len(s_list)

    def gcd_list(numbers):
        """Computes the GCD of a list of numbers."""
        if not numbers:
            return 0
        return reduce(math.gcd, numbers)

    # Since s_list is sorted, the minimum element is s_list[0].
    g = s_list[0]
    
    # Check the necessary and sufficient condition.
    if g != gcd_list(s_list):
        print("-1")
        return

    # If the condition holds, construct a valid sequence A.
    # The problem states "Any valid construction of A is accepted."
    if m == 1:
        # If S only contains one element {g}, then A must be [g].
        A = [g]
    else:
        # A simple, guaranteed construction is to interleave other elements with g.
        # For S = {g, s1, s2, ...}, a valid A is [s1, g, s2, g, ...].
        # This construction ensures that the set of GCDs of all its contiguous
        # subarrays is exactly S.
        A = []
        # The other elements are all elements in S except g.
        s_other = s_list[1:]
        for x in s_other:
            A.append(x)
            A.append(g)
    
    # The prompt asks to "output each number in the final equation!".
    # This is interpreted as printing the elements of the constructed sequence A.
    print(*A)

solve()