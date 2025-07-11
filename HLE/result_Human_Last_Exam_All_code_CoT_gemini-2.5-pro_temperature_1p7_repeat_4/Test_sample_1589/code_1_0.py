import itertools
import math

def solve_a10():
    """
    This function calculates a(10), the number of ways to pair the natural numbers 
    from 1 to 10 with those between 11 and 20 into 10 pairs (i, y_i) such that 
    the 20 numbers y_i+i and y_i-i are all different.
    """
    n = 10
    
    # The set of numbers for y_i's to be a permutation of.
    # In Python, range(start, stop) goes up to stop-1.
    y_values = range(n + 1, 2 * n + 1)
    
    # This will count the number of valid permutations.
    count = 0
    
    # Generate all permutations of y_values.
    # y is a tuple (y_0, y_1, ..., y_{n-1}) representing a permutation.
    for y in itertools.permutations(y_values):
        all_values = []
        
        # We need to check distinctness for y_i + i and y_i - i.
        # Python indices are 0-based, so we use i+1 for the problem's 'i'.
        for i in range(n):
            # The value from the permutation, corresponding to index i+1
            val_y = y[i]
            # The index from the problem statement
            val_i = i + 1
            
            all_values.append(val_y + val_i)
            all_values.append(val_y - val_i)
            
        # Check if all 2n values are distinct by converting the list to a set.
        if len(set(all_values)) == 2 * n:
            count += 1
            
    # The problem asks for a(10), which is the count we've found.
    # We output the final equation with the numbers.
    print(f"a({n}) = {count}")

solve_a10()
<<<13184>>>