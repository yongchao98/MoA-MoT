import itertools

def find_fibonacci_sum_combinations():
    """
    This function identifies all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number.

    The constraints are:
    - The combination must consist of 3, 4, 5, 6, or 7 Fibonacci numbers.
    - Fibonacci numbers cannot be duplicated, except for '1' which can be used twice.
    - The maximum Fibonacci number considered is the 40th term.
    """
    
    # Step 1: Generate Fibonacci numbers up to the 40th term.
    # The pool starts with two '1's to allow for its dual use.
    n = 40
    fib_pool = [1, 1]
    a, b = 1, 1
    for _ in range(n - 2):
        a, b = b, a + b
        fib_pool.append(b)

    # Step 2: Create a set for efficient lookup of Fibonacci numbers.
    fib_lookup_set = set(fib_pool)

    # Step 3 & 5: Iterate through specified combination sizes and store unique valid results.
    combination_sizes = range(3, 8)  # Sizes 3, 4, 5, 6, 7
    valid_combinations_set = set()

    for size in combination_sizes:
        # Step 4: Find all combinations of the current size from the Fibonacci pool.
        for combo in itertools.combinations(fib_pool, size):
            # Calculate the sum and check if it's a Fibonacci number.
            current_sum = sum(combo)
            if current_sum in fib_lookup_set:
                # Add the tuple to a set to automatically handle uniqueness.
                # itertools.combinations produces tuples where elements are sorted
                # based on their position in the input iterable.
                valid_combinations_set.add(combo)

    # Step 6: Print the results in a clear format.
    # Sort the results for a consistent and readable output, first by the sum, then by the combination itself.
    sorted_combinations = sorted(list(valid_combinations_set), key=lambda c: (sum(c), c))

    print("Found the following combinations:\n")
    for combo in sorted_combinations:
        # Build and print the equation string for each valid combination.
        equation_str = " + ".join(map(str, combo))
        total_sum = sum(combo)
        print(f"{equation_str} = {total_sum}")

    # Print the total count of unique combinations found.
    count = len(valid_combinations_set)
    print(f"\nTotal number of unique combinations found: {count}")


if __name__ == '__main__':
    find_fibonacci_sum_combinations()