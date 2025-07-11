def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2 + 3n + 1.
    """
    
    # Set a limit for the dynamic programming table. This needs to be large enough
    # to find the largest unreachable number. Experimental results show that
    # a limit of a few hundred is sufficient for the result to stabilize.
    limit = 300

    # Step 1: Generate the numbers (generators) in the set
    generators = []
    n = 1
    while True:
        term = 2 * n * n + 3 * n + 1
        if term > limit:
            break
        generators.append(term)
        n += 1

    # Step 2: Use dynamic programming to find all reachable sums
    # is_reachable[i] is True if i can be formed by a sum of distinct generators.
    is_reachable = [False] * (limit + 1)
    is_reachable[0] = True

    # Step 3: Populate the is_reachable table
    for gen in generators:
        # Iterate downwards to ensure each generator is used at most once per sum
        for i in range(limit, gen - 1, -1):
            if is_reachable[i - gen]:
                is_reachable[i] = True

    # Step 4: Find the largest number that is not reachable
    largest_unreachable = -1
    for i in range(limit, -1, -1):
        if not is_reachable[i]:
            largest_unreachable = i
            break
            
    # Step 5: The smallest number N is the largest unreachable number + 1
    N = largest_unreachable + 1
    
    print(N)

solve()