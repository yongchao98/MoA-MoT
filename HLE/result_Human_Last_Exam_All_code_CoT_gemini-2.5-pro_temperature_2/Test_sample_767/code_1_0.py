import collections

def solve():
    """
    Calculates the exact value of the limit by summing the contributions from all
    coefficient tuples that lead to a number of solutions growing linearly with N.
    """

    # Case 1: a=b=c=d=e=f=0. Equation is F_n = -g.
    # We count the number of values of n for each possible value of -g.
    
    # Generate Fibonacci numbers and count solutions for n for a given value
    # F_0=0, F_1=1, F_2=1, F_3=2, F_4=3, F_5=5, F_6=8, F_7=13, F_8=21
    # For a value V, n_sols_map[V] is the number of n such that F_n=V
    n_sols_map = collections.defaultdict(int)
    n_sols_map[0] = 1 # n=0
    n_sols_map[1] = 2 # n=1, 2
    n_sols_map[2] = 1 # n=3
    n_sols_map[3] = 1 # n=4
    n_sols_map[5] = 1 # n=5
    n_sols_map[8] = 1 # n=6
    n_sols_map[13] = 1 # n=7
    n_sols_map[21] = 1 # n=8
    
    limit_from_case1 = 0
    print("Analysis for Case 1 (a=b=c=d=e=f=0):")
    print("-" * 30)

    # Loop over g in [-25, 25]
    for g in range(-25, 26):
        k = -g
        if k in n_sols_map:
            num_n_solutions = n_sols_map[k]
            limit_from_case1 += num_n_solutions
            # For each such n, any m in [0, N-1] is a solution.
            # Number of solutions = num_n_solutions * N. Limit contribution = num_n_solutions.
            print(f"For g={g:3}, F_n = {k:2}. This provides {num_n_solutions} solution(s) for n.")
            print(f"    Contribution to the limit: {num_n_solutions}")

    print(f"\nTotal contribution from Case 1: {limit_from_case1}")
    
    # Case 2: F_n = F_m
    # This comes from coefficients a=b=c=d=e=g=0 and f=-1.
    # Number of solutions is N (for n=m) + 2 (for (1,2) and (2,1)).
    # Limit of (N+2)/N is 1.
    limit_from_case2 = 1
    
    print("\nAnalysis for Case 2 (a=b=c=d=e=0, f=-1, g=0):")
    print("-" * 30)
    print("Equation is F_n = F_m.")
    print("This gives N+2 solutions for large N.")
    print(f"Contribution to the limit: {limit_from_case2}")

    # Total limit is the sum of contributions.
    total_limit = limit_from_case1 + limit_from_case2
    print("\n" + "=" * 30)
    print(f"The final equation for the limit is the sum of all contributions:")
    print(f"Limit = (Contribution from Case 1) + (Contribution from Case 2)")
    print(f"Limit = {limit_from_case1} + {limit_from_case2} = {total_limit}")
    print("=" * 30)
    
solve()
