import math
from functools import reduce

def solve():
    """
    Solves the problem of reconstructing sequence A from set S.
    """
    
    # It's not specified how the input is given, but competitive programming
    # problems usually give the size of the set/list first.
    # Let's assume the first input is the number of elements in S,
    # and the second line is the space-separated elements of S.
    try:
        m_str = input()
        if not m_str: # Handle empty input line
            print("-1")
            return
        m = int(m_str)
        elements_str = input()
        if not elements_str: # Handle empty input line
            print("-1")
            return
        s_list = list(map(int, elements_str.split()))
    except (ValueError, IndexError):
        # Fallback for empty or invalid input
        print("-1")
        return

    if m != len(s_list):
        # Data inconsistency
        print("-1")
        return

    # The input S is a set, so we should only consider unique elements.
    s_unique = sorted(list(set(s_list)))
    
    def gcd_list(numbers):
        """Calculates the GCD of a list of numbers."""
        if not numbers:
            return 0
        return reduce(math.gcd, numbers)

    # Step 1: Find the minimum element of S.
    min_s = s_unique[0]

    # Step 2: Check the necessary and sufficient condition.
    # The minimum element of S must be the GCD of all elements in S.
    gcd_s = gcd_list(s_unique)

    if min_s != gcd_s:
        # If the condition fails, no solution exists.
        print("-1")
    else:
        # If the condition holds, a solution is guaranteed to exist.
        # We can construct one by interleaving the minimum element with the others.
        g = min_s
        
        # Handle the edge case where S has only one element.
        if len(s_unique) == 1:
            print(1)
            print(g)
            return

        # Construction: A = [s_2, g, s_3, g, ..., s_k, g]
        # where s_i are the elements of S other than g.
        reconstructed_A = []
        for x in s_unique:
            # We don't need to add g itself in this loop,
            # as it will be interleaved.
            if x != g:
                reconstructed_A.append(x)
                reconstructed_A.append(g)

        # Output the length of the constructed sequence A.
        print(len(reconstructed_A))
        
        # Output the elements of A, separated by spaces.
        # The prompt "output each number in the final equation" is interpreted
        # as printing the elements of the constructed sequence.
        print(*reconstructed_A)

solve()
<<<J>>>