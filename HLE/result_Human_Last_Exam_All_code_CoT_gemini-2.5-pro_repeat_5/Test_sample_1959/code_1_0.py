import itertools

def find_fibonacci_sum_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers whose sum is
    also a Fibonacci number, up to the 40th Fibonacci number.
    """
    
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The problem specifies the sequence starts with 1, 1.
    fib_numbers_pool = [1, 1]
    a, b = 1, 1
    # We already have the first 2, so we generate the next 38.
    for _ in range(38):
        a, b = b, a + b
        fib_numbers_pool.append(b)

    # Use a set for efficient checking of whether a sum is a Fibonacci number.
    fib_lookup_set = set(fib_numbers_pool)
    
    found_combinations = []
    
    # Step 2: Iterate through required combination lengths (3 to 7).
    for k in range(3, 8):
        # Step 3: Generate all unique combinations of size k.
        # `itertools.combinations` handles the uniqueness constraint.
        # Since fib_numbers_pool contains two '1's, combinations can include up to two '1's.
        for combo in itertools.combinations(fib_numbers_pool, k):
            
            # Step 4: Calculate the sum of the current combination.
            current_sum = sum(combo)
            
            # Step 5: Check if the sum is also a Fibonacci number.
            if current_sum in fib_lookup_set:
                # Step 6: Store the valid combination as a formatted string.
                equation = f"{' + '.join(map(str, combo))} = {current_sum}"
                found_combinations.append(equation)

    # Step 7: Print all the found equations.
    for eq in found_combinations:
        print(eq)
        
    # Print the total count of valid combinations.
    print(f"\nFound {len(found_combinations)} such combinations in total.")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()