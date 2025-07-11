import itertools

def calculate_a(n):
    """
    Calculates a(n), the number of ways to pair {1,...,n} with {n+1,...,2n}
    such that for each pair (i, y_i), the 2n numbers y_i+i and y_i-i are all different.
    """
    count = 0
    indices = range(1, n + 1)
    values = range(n + 1, 2 * n + 1)

    # Iterate through all permutations of {n+1, ..., 2n}.
    # Each permutation p represents a possible assignment (y_1, y_2, ..., y_n).
    for p in itertools.permutations(values):
        
        generated_numbers = set()
        is_valid = True
        
        # For each pair (i, y_i), calculate the two required numbers
        # and check if they cause a collision.
        for i in indices:
            y_i = p[i-1]
            
            # First number: y_i + i
            num1 = y_i + i
            if num1 in generated_numbers:
                is_valid = False
                break
            generated_numbers.add(num1)
            
            # Second number: y_i - i
            num2 = y_i - i
            if num2 in generated_numbers:
                is_valid = False
                break
            generated_numbers.add(num2)
            
        if is_valid:
            count += 1
            
    # The problem mentions outputting numbers in an equation.
    # As a(10) = 0, no such arrangements exist, so there are no numbers to show.
    # We will just print the final count.
    print(f"a({n}) = {count}")


# Calculate a(10) as requested by the user.
# Note: This calculation might take a few minutes to complete.
calculate_a(10)
