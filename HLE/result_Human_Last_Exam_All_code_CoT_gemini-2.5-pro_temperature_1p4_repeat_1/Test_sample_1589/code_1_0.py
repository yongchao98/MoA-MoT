import itertools

def solve_a(n):
    """
    Calculates a(n), the number of ways to pair {1..n} with {n+1..2n}
    such that for each pair (x_i, y_i), the 2n values y_i+x_i and y_i-x_i
    are all distinct.
    """
    # The set {x_1, ..., x_n} is {1, ..., n}. We can fix x_i = i.
    # The set {y_1, ..., y_n} is a permutation of {n+1, ..., 2n}.
    # A permutation p of {1, ..., n} can define the pairing: y_i = n + p_i.
    
    count = 0
    # Generate all permutations of the numbers (1, 2, ..., n)
    # The permutation p is a tuple where p[i-1] corresponds to p(i).
    for p in itertools.permutations(range(1, n + 1)):
        # This set will store the 2n values to check for distinctness.
        all_values = set()
        is_valid_permutation = True
        
        # For each i from 1 to n, calculate the two values.
        for i in range(1, n + 1):
            x_i = i
            p_i = p[i-1]  # p(i)
            y_i = n + p_i
            
            val_plus = y_i + x_i
            val_minus = y_i - x_i
            
            # Check for duplicates. If a duplicate is found, this permutation is invalid.
            if val_plus in all_values:
                is_valid_permutation = False
                break
            all_values.add(val_plus)
            
            if val_minus in all_values:
                is_valid_permutation = False
                break
            all_values.add(val_minus)
            
        # If the inner loop completed without a break, all 2n values were distinct.
        if is_valid_permutation:
            count += 1
            
    return count

if __name__ == '__main__':
    n = 10
    result = solve_a(n)
    # The final equation is a(10) = result
    print(f"a({n}) = {result}")
