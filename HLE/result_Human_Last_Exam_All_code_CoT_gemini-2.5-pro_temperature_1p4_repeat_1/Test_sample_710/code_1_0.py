def find_frobenius_number_for_distinct_sum():
    """
    This function determines the smallest number N such that any number M >= N
    can be written as a sum of distinct numbers of the form f(n) = 2n^2 + 3n + 1.

    As explained in the plan, the set of generators f(n) for all integers n is
    the set of all triangular numbers {1, 3, 6, 10, ...}.
    This function finds the largest number that cannot be written as a sum of
    distinct triangular numbers. N is that number plus one.
    """
    
    # We check for reachable sums up to a reasonable limit.
    # Literature suggests the largest unreachable number is small, so 100 is sufficient.
    MAX_SUM = 100

    # Generate triangular numbers (our set of generators) up to MAX_SUM.
    generators = []
    k = 1
    while True:
        t_k = k * (k + 1) // 2
        if t_k > MAX_SUM:
            break
        generators.append(t_k)
        k += 1

    # reachable[i] will be true if i can be formed by a sum of distinct generators.
    # We initialize reachable[0] as True for the empty sum.
    reachable = [False] * (MAX_SUM + 1)
    reachable[0] = True

    # Use dynamic programming to solve the subset sum problem.
    for gen in generators:
        # Iterate backwards to ensure each generator is used at most once for any given sum.
        for i in range(MAX_SUM, gen - 1, -1):
            if reachable[i - gen]:
                reachable[i] = True

    # Find the largest number up to MAX_SUM that is not reachable.
    largest_unreachable = -1
    for i in range(MAX_SUM, -1, -1):
        if not reachable[i]:
            largest_unreachable = i
            break
    
    # N is the smallest number such that all numbers >= N are representable.
    # This is equal to the largest unreachable number + 1.
    N = largest_unreachable + 1
    
    print(f"The largest number that cannot be written as a sum of distinct numbers of the form 2n^2+3n+1 is {largest_unreachable}.")
    print("The smallest number N such that any number >= N can be written as such a sum is found by the equation:")
    # Here we print the numbers in the final equation that determines N.
    print(f"{N} = {largest_unreachable} + 1")

find_frobenius_number_for_distinct_sum()
<<<34>>>