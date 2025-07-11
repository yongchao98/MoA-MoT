import math

def combinations(n, k):
    """Calculates the number of combinations C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_main_diagonal_symmetry(n):
    """
    Calculates the number of configurations symmetric along the main diagonal.
    This is the number of involutions a_n.
    Recurrence: a_n = a_{n-1} + (n-1) * a_{n-2}
    """
    a = [0] * (n + 1)
    a[0] = 1
    if n >= 1:
        a[1] = 1
    for i in range(2, n + 1):
        a[i] = a[i - 1] + (i - 1) * a[i - 2]
    return a[n]

def calculate_anti_diagonal_symmetry(n):
    """
    Calculates the number of configurations symmetric along the anti-diagonal (for even n).
    This is the number b_n.
    Recurrence: b_n = 2*b_{n-2} + (n-2)*b_{n-4}
    """
    if n % 2 != 0:
        return 0
    b = [0] * (n + 1)
    b[0] = 1
    if n >= 2:
        b[2] = 2
    for i in range(4, n + 1, 2):
        b[i] = 2 * b[i - 2] + (i - 2) * b[i - 4]
    return b[n]

def calculate_both_diagonals_symmetry(n):
    """
    Calculates the number of configurations symmetric along both diagonals.
    This is the number of involutions on n/2 sets X_k = {k, n+1-k},
    summing 2^c(P) over all such involutions P, where c(P) is the number of cycles in P.
    """
    if n % 2 != 0:
        return 0
    k = n // 2
    total = 0
    # j is the number of 2-cycles in the involution P on k sets
    for j in range(k // 2 + 1):
        # Number of ways to choose an involution on k items with j 2-cycles
        # This is C(k, 2j) * (2j-1)!!
        try:
            num_involutions = combinations(k, 2 * j)
            # (2j-1)!! = (2j)! / (j! * 2^j)
            if j > 0:
              num_involutions *= math.factorial(2 * j) // (math.factorial(j) * (2 ** j))
        except (ValueError, IndexError):
            num_involutions = 0

        # Number of cycles in P is (k-2j) fixed points + j 2-cycles = k-j
        num_cycles = k - j
        total += num_involutions * (2 ** num_cycles)
    return total


def solve():
    """
    Solves the checkerboard problem for an 8x8 board.
    """
    n = 8
    
    # Calculate N(main): number of configurations symmetric along the main diagonal
    num_main = calculate_main_diagonal_symmetry(n)
    print(f"Number of placements symmetric along the main diagonal: {num_main}")

    # Calculate N(anti): number of configurations symmetric along the anti-diagonal
    num_anti = calculate_anti_diagonal_symmetry(n)
    print(f"Number of placements symmetric along the anti-diagonal: {num_anti}")

    # Calculate N(both): number of configurations symmetric along both diagonals
    num_both = calculate_both_diagonals_symmetry(n)
    print(f"Number of placements symmetric along both diagonals: {num_both}")

    # Use the Principle of Inclusion-Exclusion to find the total number of configurations
    total_configs = num_main + num_anti - num_both
    
    print("\nThe problem asks for the number of configurations symmetric along at least one diagonal.")
    print("Using the Principle of Inclusion-Exclusion:")
    print(f"Total = N(main) + N(anti) - N(both)")
    print(f"Total = {num_main} + {num_anti} - {num_both} = {total_configs}")

solve()