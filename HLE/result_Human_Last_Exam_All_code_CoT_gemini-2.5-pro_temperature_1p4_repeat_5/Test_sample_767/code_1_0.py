import collections

def solve_limit():
    """
    This script calculates the exact value of the limit based on the analytical solution.
    The limit is equal to the number of solutions (n, g) for the equation F_n = -g,
    under the problem's constraints. This is equivalent to counting the number of solutions
    n for F_n = V, where V is an integer from 0 to 25.
    """
    
    V_max = 25
    
    # Create a map from a Fibonacci value to the list of indices 'n' that produce it.
    fib_map = collections.defaultdict(list)
    
    # Handle F_0 and F_1
    fib_map[0].append(0)
    fib_map[1].append(1)
    
    f_prev, f_curr = 0, 1
    n = 2
    while True:
        f_next = f_prev + f_curr
        if f_next > V_max:
            break
        fib_map[f_next].append(n)
        f_prev, f_curr = f_curr, f_next
        n += 1

    total_solutions_K = 0
    solution_counts = []
    
    # We are interested in any Fibonacci number that falls within the range [0, 25].
    # These are the only values of V = -g that yield solutions for n.
    # The Fibonacci numbers in this range are 0, 1, 2, 3, 5, 8, 13, 21.
    
    fib_keys_in_range = sorted([k for k in fib_map.keys() if k <= V_max])
    
    for v_val in fib_keys_in_range:
        count = len(fib_map[v_val])
        total_solutions_K += count
        solution_counts.append(str(count))

    # The problem asks to output each number in the final equation.
    # The final "equation" is the sum of the counts of solutions.
    final_equation_str = " + ".join(solution_counts)

    print("The limit is the sum of the number of solutions for n in F_n = V for each Fibonacci number V in [0, 25]:")
    print(f"{final_equation_str} = {total_solutions_K}")

solve_limit()