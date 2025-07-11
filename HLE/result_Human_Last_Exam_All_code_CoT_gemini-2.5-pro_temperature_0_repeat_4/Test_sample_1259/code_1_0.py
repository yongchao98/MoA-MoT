import math

def solve():
    """
    Reads a set S and attempts to restore a sequence A.
    """
    try:
        # In some environments, input() might include the size of the set on the first line.
        # This code assumes the first line of input is the list of numbers.
        # If the first line is the count 'm', it should be read and ignored.
        # For example:
        # m = int(input())
        # s_list = list(map(int, input().split()))
        
        s_list = list(map(int, input().split()))
    except (IOError, ValueError):
        # Handle potential empty input or malformed input
        print("-1")
        return

    if not s_list:
        print("-1")
        return

    # The input is a set, so we should only consider unique elements
    s_set = set(s_list)
    s_list = sorted(list(s_set))
    
    m = len(s_list)
    g = s_list[0]

    # Check the necessary and sufficient condition
    for x in s_list:
        if x % g != 0:
            print("-1")
            return

    # If the condition holds, construct a sequence A
    if m == 1:
        # If S has only one element, A is just that element
        print(1)
        print(g)
    else:
        # Use the interleaving construction for |S| > 1
        # A = [s_2, g, s_3, g, ..., s_m, g]
        result_a = []
        # s_list is sorted, so s_list[0] is g.
        # We interleave the other elements with g.
        for i in range(1, m):
            result_a.append(s_list[i])
            result_a.append(g)
        
        print(len(result_a))
        # The problem asks to output each number in the final sequence
        print(*result_a)

solve()