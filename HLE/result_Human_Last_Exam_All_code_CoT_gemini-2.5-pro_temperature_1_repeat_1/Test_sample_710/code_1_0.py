import sys

def solve():
    """
    Determines the smallest number N such that any number >= N can be written
    as a sum of distinct numbers of the form 2n^2+3n+1.
    """
    # The upper limit for our search. The numbers in the set grow quadratically,
    # and the sum of the first k numbers grows cubically. This suggests that
    # the largest unreachable number will not be excessively large. A limit of
    # 400 is sufficient.
    limit = 400

    # Generate the numbers s_n = (n+1)(2n+1)
    generators = []
    n = 1
    while True:
        s = (n + 1) * (2 * n + 1)
        if s > limit:
            break
        generators.append(s)
        n += 1

    # Use dynamic programming to find all reachable sums
    is_reachable = [False] * (limit + 1)
    is_reachable[0] = True
    
    # Store predecessors to reconstruct the sum
    predecessor = {}

    for s in generators:
        for j in range(limit, s - 1, -1):
            if is_reachable[j - s] and not is_reachable[j]:
                is_reachable[j] = True
                predecessor[j] = s

    # Find the largest unreachable number
    last_unreachable = 0
    for i in range(1, limit + 1):
        if not is_reachable[i]:
            last_unreachable = i
    
    n_value = last_unreachable + 1
    print(f"The smallest number N is: {n_value}")
    print("-" * 30)
    print(f"Any number >= {n_value} can be written as a sum of distinct numbers of the form 2n^2+3n+1.")
    print("Here are some examples:")

    def get_decomposition(num):
        """Finds the sum decomposition for a given number."""
        if num not in predecessor and num != 0:
            return None
        parts = []
        curr = num
        while curr > 0:
            part = predecessor[curr]
            parts.append(part)
            curr -= part
        return sorted(parts)

    # Print the decomposition for N, N+1, N+2, and N+3
    for i in range(n_value, n_value + 4):
        decomposition = get_decomposition(i)
        if decomposition:
            equation = f"{i} = {' + '.join(map(str, decomposition))}"
            print(equation)
        else:
            print(f"Could not find decomposition for {i}")

solve()