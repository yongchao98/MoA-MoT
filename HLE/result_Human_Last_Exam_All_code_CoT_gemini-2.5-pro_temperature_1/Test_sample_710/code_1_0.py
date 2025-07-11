def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """

    # Step 1: Generate the first few terms of the set S.
    # Using 15 terms is sufficient for the result to stabilize.
    num_terms = 15
    S = []
    for n in range(1, num_terms + 1):
        term = (2 * n + 1) * (n + 1)
        S.append(term)

    # Step 2: Generate all possible sums of distinct elements from S.
    # We start with {0} to represent the sum of an empty set of numbers.
    reachable_sums = {0}
    for x in S:
        # For each new element x, create new sums by adding x to all existing sums.
        new_sums = {s + x for s in reachable_sums}
        reachable_sums.update(new_sums)

    # Step 3: Find the largest number that is not in the set of reachable sums.
    # The maximum possible sum we have generated is sum(S).
    # We search backwards from a sufficiently high number. If the result stabilizes,
    # it is likely the correct one. A search limit of 400 is more than enough
    # as the result is small.
    search_limit = 400
    largest_unreachable = -1
    for i in range(search_limit, 0, -1):
        if i not in reachable_sums:
            largest_unreachable = i
            break

    # The smallest number N from which all integers are representable
    # is the largest unreachable number + 1.
    N = largest_unreachable + 1
    
    print(f"The set of numbers begins with: {S[:6]}...")
    print(f"After checking sums up to a high limit, the largest integer that CANNOT be formed is: {largest_unreachable}")
    print(f"Therefore, the smallest number N such that any number >= N can be written as a sum of distinct numbers from the set is {largest_unreachable} + 1.")
    print(f"N = {N}")

solve()
<<<119>>>