import math

def solve():
    """
    Solves the problem by first identifying the building blocks as triangular numbers,
    then using dynamic programming to find the largest number not representable as their distinct sum.
    """
    
    # Step 1 & 2: Identify numbers and restate the problem.
    # The numbers of the form 2n^2+3n+1 for integer n are precisely the set of
    # triangular numbers T_k = k(k+1)/2 for k >= 1.
    # The problem is to find the largest number M not representable as a sum
    # of distinct triangular numbers. The answer N will be M + 1.

    # Step 3: Use dynamic programming to find all reachable sums up to a limit.
    # A limit of 100 should be sufficient to find the pattern of unreachability.
    limit = 100
    is_reachable = [False] * (limit + 1)
    is_reachable[0] = True
    
    triangulars = []
    k = 1
    while True:
        t = k * (k + 1) // 2
        if t > limit:
            break
        triangulars.append(t)
        # For each new triangular number, update the reachable sums
        for s in range(limit, t - 1, -1):
            if is_reachable[s - t]:
                is_reachable[s] = True
        k += 1

    # Find the largest unreachable number M
    unreachable_numbers = [i for i, reachable in enumerate(is_reachable) if not reachable]
    
    if not unreachable_numbers:
        print("Error: No unreachable numbers found. The limit might be too low.")
        return

    largest_unreachable = unreachable_numbers[-1]
    N = largest_unreachable + 1
    
    print(f"The numbers of the form 2n^2+3n+1 are the triangular numbers: {triangulars[:7]}...")
    print(f"The numbers that cannot be written as a sum of distinct triangular numbers are:")
    print(", ".join(map(str, unreachable_numbers)))
    print(f"\nThe largest number not representable is {largest_unreachable}.")
    print(f"Therefore, the smallest number N such that any number >= N is representable is {N}.")

    # Step 4: Decompose N to satisfy the output format requirement.
    # We find the representation of N = 34. A greedy approach works here.
    # 34 = 28 + 6.
    # 28 is T_7. For odd k=7, n = (k-1)/2 = 3. So 28 = f(3).
    # 6 is T_3. For odd k=3, n = (k-1)/2 = 1. So 6 = f(1).
    
    print("\nTo demonstrate, we show the decomposition for N = 34:")
    n1, n2 = 3, 1
    val1 = 2 * n1**2 + 3 * n1 + 1
    val2 = 2 * n2**2 + 3 * n2 + 1
    
    print(f"{N} = {val1} + {val2}")
    print(f"{N} = (2*{n1}^2 + 3*{n1} + 1) + (2*{n2}^2 + 3*{n2} + 1)")


solve()
<<<34>>>