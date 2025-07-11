import sys

# Set a higher recursion limit for deep calculations, although this iterative approach doesn't need it.
# sys.setrecursionlimit(20000)

def solve_limit_ratio():
    """
    Calculates the limit of the expected ratio of remaining numbers to n.
    """
    # N is the large number up to which we calculate the expectation.
    N = 10000

    # E is a list to store the expected values E_n.
    E = [0.0] * (N + 1)

    # Base cases from the problem definition.
    # E_0 = 0 (no numbers)
    # E_1 = 1 (one number, which is isolated)
    E[0] = 0
    if N >= 1:
        E[1] = 1
    
    # E_2 can be calculated from the sum formula: E_2 = 2/(2-1) * E_0 = 0.
    if N >= 2:
        E[2] = 0

    # Iteratively compute E_n using the recurrence relation
    # n*E_{n+1} = (n-1)*E_n + 2*E_{n-1} which we rewrite for E_n
    # as (n-1)*E_n = (n-2)*E_{n-1} + 2*E_{n-2} for n >= 3.
    for n in range(3, N + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # The ratio we want to find the limit of.
    ratio = E[N] / N

    # As requested, output the final equation used for the last step of the calculation.
    print(f"The recurrence relation for E_n is: (n-1)*E_n = (n-2)*E_{{n-1}} + 2*E_{{n-2}}")
    print(f"Calculation for the final term E[{N}]:")
    print(f"({N}-1)*E[{N}] = ({N}-2)*E[{N}-1] + 2*E[{N}-2]")
    print(f"{N-1}*E[{N}] = {N-2}*{E[N-1]:.4f} + 2*{E[N-2]:.4f}")
    val1 = (N-2) * E[N-1]
    val2 = 2 * E[N-2]
    print(f"{N-1}*E[{N}] = {val1:.4f} + {val2:.4f}")
    sum_val = val1 + val2
    print(f"{N-1}*E[{N}] = {sum_val:.4f}")
    final_E_N = sum_val / (N - 1)
    print(f"E[{N}] = {final_E_N:.4f}")
    
    print("\nLimit of the ratio E_n / n as n -> infinity:")
    print(f"For n = {N}, the ratio E_n / n is {ratio}")
    print(f"The limit appears to converge to 1/6.")
    
solve_limit_ratio()