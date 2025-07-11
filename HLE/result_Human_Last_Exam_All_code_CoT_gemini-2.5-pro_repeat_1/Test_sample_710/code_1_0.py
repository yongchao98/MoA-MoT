def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """
    # We need to find the largest number that cannot be formed. N will be that number + 1.
    # We use dynamic programming to find all reachable sums up to a certain limit.
    LIMIT = 200  # A sufficiently large number to find the turning point.

    # Generate the set of numbers of the form 2n^2+3n+1.
    # The formula can be factored to (2n+1)(n+1).
    nums = []
    n = 0
    while True:
        val = (2 * n + 1) * (n + 1)
        if val >= LIMIT:
            break
        # We only consider positive numbers from non-negative n.
        if val > 0:
            nums.append(val)
        n += 1

    # reachable[i] will be true if i can be expressed as a sum of distinct nums.
    reachable = [False] * LIMIT
    reachable[0] = True  # Sum of 0 is possible with an empty set.

    # For each number, update the reachable sums.
    # Iterate backwards to ensure each number is used at most once (distinct sums).
    for num in nums:
        for i in range(LIMIT - 1, num - 1, -1):
            if reachable[i - num]:
                reachable[i] = True

    # Find the largest number that is not reachable.
    largest_unreachable = -1
    for i in range(LIMIT - 1, 0, -1):
        if not reachable[i]:
            largest_unreachable = i
            break

    # N is the first number such that all numbers >= N are reachable.
    # This is largest_unreachable + 1.
    if largest_unreachable != -1:
        N = largest_unreachable + 1
        print(f"The largest unreachable number is: {largest_unreachable}")
        print(f"The smallest number N is the next integer.")
        print(f"The final equation is: {largest_unreachable} + 1 = {N}")
    else:
        print(f"All numbers up to {LIMIT} seem to be reachable.")

solve()