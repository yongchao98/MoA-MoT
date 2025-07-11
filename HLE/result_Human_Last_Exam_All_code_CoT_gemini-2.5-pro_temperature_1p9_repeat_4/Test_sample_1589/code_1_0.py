import itertools

def calculate_a(n):
    """
    Calculates a(n), the number of ways to pair numbers {1..n} with {n+1..2n}
    such that for each pair (i, y_i), the 2n values y_i+i and y_i-i are all different.
    """
    
    # The set of numbers to be permuted for the y_i values.
    y_domain = range(n + 1, 2 * n + 1)
    
    # Initialize a counter for the number of valid ways.
    count = 0
    
    # Generate all permutations of y_domain. Each permutation is a tuple
    # representing a possible assignment (y_1, y_2, ..., y_n).
    for y_permutation in itertools.permutations(y_domain):
        
        generated_values = set()
        is_valid_permutation = True
        
        # Iterate through i from 1 to n.
        for i in range(1, n + 1):
            # Get the corresponding y_i for the current i (note: y_permutation is 0-indexed).
            y_i = y_permutation[i-1]
            
            # Calculate the two values for the current pair.
            val_plus = y_i + i
            val_minus = y_i - i
            
            # Check if either value causes a collision with previously generated values.
            if val_plus in generated_values or val_minus in generated_values:
                is_valid_permutation = False
                break
            
            # If no collision, add the new values to the set.
            generated_values.add(val_plus)
            generated_values.add(val_minus)
            
        # If the inner loop completed without a break, the permutation is valid.
        if is_valid_permutation:
            count += 1
            
    return count

def main():
    """
    Main function to calculate and print a(10).
    """
    n = 10
    result = calculate_a(n)
    
    # As requested, output the number in a final equation format.
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    main()
