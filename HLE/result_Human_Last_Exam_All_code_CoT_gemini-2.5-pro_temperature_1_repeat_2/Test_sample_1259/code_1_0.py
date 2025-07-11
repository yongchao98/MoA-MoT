import math
import sys

def solve():
    """
    Reads a set S, and if possible, restores a sequence A that generates S.
    The logic is based on the necessary and sufficient condition that min(S)
    must divide all elements of S.
    """
    try:
        # Read the single line of space-separated integers for the set S
        line = sys.stdin.readline().strip()
        if not line:
            return
        s_list = list(map(int, line.split()))
    except (IOError, ValueError):
        # Handle potential read errors or non-integer input
        print("-1")
        return

    if not s_list:
        print("-1")
        return
    
    # Use a set for efficient lookup and to handle potential duplicates in input
    s_set = set(s_list)
    
    # Find the minimum element of the set S
    g = min(s_set)
    
    # Check the necessary condition: min(S) must divide all elements of S.
    # If this condition fails, no solution exists.
    for x in s_set:
        if x % g != 0:
            print("-1")
            return
            
    # The condition holds, which is also sufficient. A solution is guaranteed.
    # We can now construct a valid sequence A. A provably correct (though
    # not always shortest) construction is to interleave g with the other
    # elements of S.
    
    # Get the other elements of S, sorted for consistent output
    s_other = sorted([x for x in s_set if x != g])
    
    # Construct the sequence A = [g, s_2, g, s_3, g, ..., s_m, g]
    result_a = [g]
    for val in s_other:
        result_a.append(val)
        result_a.append(g)

    # Print the restored sequence A. The phrasing "output each number
    # in the final equation!" is interpreted as printing the elements of A.
    print(*result_a)

solve()
<<<J>>>