import math

def solve():
    """
    Solves the Restoring Sequence problem.
    """
    try:
        # Read input from stdin
        line1 = input().strip()
        if not line1:
            return
            
        m = int(line1)
        s_list = list(map(int, input().strip().split()))

    except (IOError, ValueError):
        # Handle cases with no input or malformed input
        # In a competitive programming context, we might assume valid input.
        print("-1")
        return

    if m == 0:
        print("0")
        print("")
        return

    # Step 1: Check the necessary and sufficient condition.
    # The minimum element of S must be the gcd of all elements in S.
    
    # Find the gcd of all elements in the set S
    if not s_list:
        # Should not happen given constraints m >= 1
        print("-1")
        return

    gcd_s = s_list[0]
    for i in range(1, m):
        gcd_s = math.gcd(gcd_s, s_list[i])

    # Find the minimum element
    min_s = min(s_list)
    
    # Check the condition
    if gcd_s != min_s:
        print("-1")
        return

    # Step 2: The condition holds, so a solution exists. Construct one.
    # We use the interleaved construction: [g, s1, g, s2, ...]
    g = min_s
    result_a = []
    
    # Add the first element, g.
    result_a.append(g)

    # Alternate other elements with g
    for x in s_list:
        if x != g:
            result_a.append(x)
            result_a.append(g)

    # In the case where S only has one element (e.g., S={5})
    if len(s_list) == 1:
        result_a = [s_list[0]]
        
    # Step 3: Print the result in the specified format
    print(len(result_a))
    # Use print with '*' to unpack the list into arguments
    print(*result_a)

solve()
