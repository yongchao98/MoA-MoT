import math
import sys

def solve():
    """
    Reads a set of integers S, checks if a generating sequence A can exist,
    and if so, constructs and prints one such A.
    """
    try:
        # Read the number of elements in S (m) and the elements themselves.
        # This part is for standard competitive programming input format.
        m_str = sys.stdin.readline()
        if not m_str: return
        m = int(m_str)
        if m == 0:
            print(0)
            print()
            return
        s_elements = list(map(int, sys.stdin.readline().split()))
    except (IOError, ValueError):
        # Fallback for local testing if input reading fails.
        print("-1")
        return

    # The problem is about a set S, so we only consider unique elements.
    s_set = sorted(list(set(s_elements)))
    
    # Find the minimum element of the set.
    min_s = s_set[0]

    # Calculate the greatest common divisor (GCD) of all elements in the set.
    if len(s_set) == 1:
        gcd_s = s_set[0]
    else:
        gcd_s = math.gcd(s_set[0], s_set[1])
        for i in range(2, len(s_set)):
            gcd_s = math.gcd(gcd_s, s_set[i])

    # A solution exists if and only if the minimum element is the GCD of the set.
    if min_s != gcd_s:
        print("-1")
    else:
        # If a solution exists, we can construct one.
        # A simple valid construction is to interleave the minimum element
        # with all other elements of the set.
        # For S = {s1, s2, s3}, A = [s1, s2, s1, s3]
        
        # Handle the case of a single-element set.
        if len(s_set) == 1:
            constructed_a = [min_s]
        else:
            constructed_a = []
            # Start with the minimum element.
            constructed_a.append(min_s)
            # Append each other element followed by the minimum element.
            for element in s_set:
                if element != min_s:
                    constructed_a.append(element)
                    constructed_a.append(min_s)
        
        # Print the length of the constructed sequence A.
        print(len(constructed_a))
        # Print the elements of the constructed sequence A.
        print(*constructed_a)

solve()