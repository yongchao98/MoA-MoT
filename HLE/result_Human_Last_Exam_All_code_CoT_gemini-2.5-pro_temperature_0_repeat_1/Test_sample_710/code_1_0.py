def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2 + 3n + 1.
    """
    # Set a sufficiently large limit for our search space.
    # A stable result for the largest unreachable number with a limit of 300
    # gives high confidence in the answer.
    MAX_LIMIT = 300

    # 1. Generate the elements of the set S up to MAX_LIMIT.
    # The formula is s_n = 2n^2 + 3n + 1 = (2n + 1)(n + 1).
    s_elements = []
    n = 0
    while True:
        val = (2 * n + 1) * (n + 1)
        if val > MAX_LIMIT:
            break
        s_elements.append(val)
        n += 1

    # 2. Use dynamic programming to find all reachable sums.
    # reachable[i] will be true if i can be formed by a sum of distinct elements.
    reachable = [False] * (MAX_LIMIT + 1)
    reachable[0] = True  # A sum of 0 is possible (by choosing no elements).

    for s in s_elements:
        # Iterate backwards to ensure each element 's' is used at most once per sum.
        for j in range(MAX_LIMIT, s - 1, -1):
            if reachable[j - s]:
                reachable[j] = True

    # 3. Find the largest number that is not reachable.
    largest_unreachable = 0
    for i in range(1, MAX_LIMIT + 1):
        if not reachable[i]:
            largest_unreachable = i

    # 4. The smallest number N is the one after the largest unreachable number.
    N = largest_unreachable + 1

    print(f"The set S begins with: {s_elements[:7]}...")
    print(f"The largest number that cannot be formed by a sum of distinct elements of S is: {largest_unreachable}")
    print(f"Therefore, the smallest number N such that any number >= N can be formed is: {N}")

solve()