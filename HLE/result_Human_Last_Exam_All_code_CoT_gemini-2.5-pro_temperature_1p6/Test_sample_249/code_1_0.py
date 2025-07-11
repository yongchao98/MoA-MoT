import math

def solve():
    """
    This function calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    I will use n=7 and m=4 as example values.
    """
    n = 7
    m = 4

    # Any tree with V > 2 vertices must have at least 2 leaves.
    # V = n+2 >= 3, so we must have m >= 2.
    # Also, the sum of degrees 2*(V-1) must be >= sum of degrees of leaves (m) + 2 * (num internal nodes).
    # 2*(n+1) >= m + 2*(n+2-m) --> 2n+2 >= m+2n+4-2m --> m >= 2.
    if m < 2 and n+2 > 2:
      print(f"A tree with {n+2} vertices cannot have {m} leaf/leaves.")
      return

    print(f"Calculating for n={n} and m={m}:")
    # Case 1: m = n + 1 (1 internal vertex)
    if m == n + 1:
        diameter = 2
        print(f"Here, m = n + 1, so there is 1 internal vertex.")
        print(f"The tree is a star graph, so the minimum diameter is 2.")
        print(f"<<<{diameter}>>>")

    # Case 2: m = n (2 internal vertices)
    elif m == n:
        diameter = 3
        print(f"Here, m = n, so there are 2 internal vertices.")
        print(f"The minimum diameter is 3.")
        print(f"<<<{diameter}>>>")
    
    # Case 3: m < n (3 or more internal vertices)
    elif m < n:
        # Subcase 3a: 2m >= n + 1
        if 2 * m >= n + 1:
            diameter = 4
            print(f"Here, m < n and 2*m >= n+1 (2*{m}={2*m} >= {n}+1={n+1}).")
            print(f"The internal vertices can form a star-like structure.")
            print(f"The minimum diameter is 4.")
            print(f"<<<{diameter}>>>")

        # Subcase 3b: 2m < n + 1
        else:
            diameter = n + 5 - 2 * m
            print(f"Here, m < n and 2*m < n+1 (2*{m}={2*m} < {n}+1={n+1}).")
            print(f"The minimum diameter is given by the formula n + 5 - 2*m.")
            print(f"Result: {n} + 5 - 2*{m} = {diameter}")
            print(f"<<<{diameter}>>>")
    else: # m > n + 1
        print("A tree with n+2 vertices and m leaves is not possible if m > n+1, as this would require <1 internal nodes.")

solve()
