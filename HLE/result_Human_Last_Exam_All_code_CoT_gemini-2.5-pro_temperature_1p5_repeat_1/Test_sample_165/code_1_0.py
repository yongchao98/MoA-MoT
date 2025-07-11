import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    # Using math.comb for efficiency and accuracy
    return math.comb(n, k)

def solve(n, k):
    """
    Calculates the number of initial sets S for which it is impossible to reach all zeros.
    n: An odd positive integer > 1.
    k: A positive integer >= n.
    """
    if not (isinstance(n, int) and n > 1 and n % 2 == 1):
        print("Error: n must be an odd positive integer greater than 1.")
        return
    if not (isinstance(k, int) and k >= n):
        print(f"Error: k must be a positive integer greater than or equal to n (k={k}, n={n}).")
        return

    # Step 1: Count odd and even numbers in [-k, k]
    # Odds in [1, k]: math.ceil(k/2)
    # Total odds: 2 * ceil(k/2) because of negative counterparts
    num_odd = 2 * math.ceil(k / 2)
    
    # Evens in [1, k]: math.floor(k/2)
    # Total evens: 2 * floor(k/2) + 1 (for 0)
    num_even = 2 * math.floor(k / 2) + 1
    
    print(f"For k = {k}:")
    print(f"Number of available odd integers (N_O) in [{-k}, {k}] = {num_odd}")
    print(f"Number of available even integers (N_E) in [{-k}, {k}] = {num_even}\n")
    
    # Step 2: Count sets with an odd number of odd integers
    # We choose j odd integers and (n-j) even integers, where j is odd.
    
    total_impossible_sets = 0
    equation_parts = []
    
    # n is odd, so we iterate through odd j = 1, 3, 5, ..., n
    for j in range(1, n + 1, 2):
        num_ways_odd = combinations(num_odd, j)
        num_ways_even = combinations(num_even, n - j)
        
        term = num_ways_odd * num_ways_even
        total_impossible_sets += term
        
        part_str = f"C({num_odd}, {j})*C({num_even}, {n-j})"
        equation_parts.append(part_str)

    print("The number of impossible sets is the sum of ways to choose j odd numbers and (n-j) even numbers, for j in {1, 3, ..., n}:")
    
    # Print the equation
    full_equation = " + ".join(equation_parts)
    print(f"Total = {full_equation}")

    # Print the values for each term in the equation
    term_values = []
    for j in range(1, n + 1, 2):
        ways_odd = combinations(num_odd, j)
        ways_even = combinations(num_even, n-j)
        term_values.append(f"{ways_odd}*{ways_even}")
    
    full_values = " + ".join(term_values)
    print(f"      = {full_values}")

    term_results = []
    for j in range(1, n + 1, 2):
        ways_odd = combinations(num_odd, j)
        ways_even = combinations(num_even, n-j)
        term_results.append(str(ways_odd * ways_even))
        
    final_sum = " + ".join(term_results)
    print(f"      = {final_sum}")
    
    print(f"      = {total_impossible_sets}")
    
    return total_impossible_sets

# Example usage with sample values, e.g., n=3, k=3
# You can change these values to solve for a specific case.
if __name__ == '__main__':
    # Per the problem description
    n = 3 
    k = 3

    solve(n, k)