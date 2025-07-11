import itertools

def solve_a_n(n):
    """
    Calculates a(n), which is the number of valid permutations for this problem.
    A permutation pi of {1, ..., n} is valid if the following conditions hold:
    1. The set F = {pi(i) + i | i in {1, ..., n}} has n distinct elements.
    2. The set G = {pi(i) - i | i in {1, ..., n}} has n distinct elements.
    3. The sets F and G are disjoint.
    """
    
    count = 0
    # We define the set of numbers to be permuted, {1, 2, ..., n}.
    numbers = range(1, n + 1)
    
    # Iterate through all permutations of the numbers {1, ..., n}.
    # p will be a tuple representing a permutation.
    # If p = (p_0, p_1, ...), then pi(i) corresponds to p[i-1].
    for p in itertools.permutations(numbers):
        
        # F corresponds to the set of {pi(i) + i} values.
        # We use p[i-1] for pi(i) and i for the index.
        f_values = {p[i-1] + i for i in numbers}
        if len(f_values) != n:
            continue
            
        # G corresponds to the set of {pi(i) - i} values.
        g_values = {p[i-1] - i for i in numbers}
        if len(g_values) != n:
            continue
            
        # Check for disjointness between F and G.
        if not f_values.isdisjoint(g_values):
            continue
            
        # If all conditions are satisfied, we've found a valid permutation.
        count += 1
        
    return count

def main():
    """
    Main function to calculate and print a(10).
    """
    n = 10
    result = solve_a_n(n)
    
    # Printing the final equation with all its numbers as requested.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    main()
