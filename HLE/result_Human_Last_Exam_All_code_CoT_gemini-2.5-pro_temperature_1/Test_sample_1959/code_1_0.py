import itertools

def find_fibonacci_sum_combinations():
    """
    This function finds and prints all combinations of 3 to 7 unique Fibonacci
    numbers whose sum is also a Fibonacci number.
    
    The constraints are:
    - Combinations must have 3, 4, 5, 6, or 7 numbers.
    - The Fibonacci numbers used must be from the first 40 terms of the sequence.
    - The number 1 can be used up to twice, as it appears twice in the sequence (F_1 and F_2).
    """

    # Step 1: Generate Fibonacci numbers.
    # We need the first 40 Fibonacci numbers (F_1 to F_40) as the pool for combinations.
    # The list will start with [1, 1, ...].
    # We also need a larger set of Fibonacci numbers to check if the sum is a Fibonacci number.
    # The maximum possible sum (F_34 + ... + F_40) is less than F_42, so generating up to
    # F_43 is sufficient for the check.
    n_max = 40
    fibs_to_generate = 43 
    
    fib_sequence = [1, 1]
    for i in range(2, fibs_to_generate):
        fib_sequence.append(fib_sequence[i-1] + fib_sequence[i-2])

    # The list of numbers to choose combinations from (F_1 to F_40).
    fib_choices = fib_sequence[:n_max]

    # A set for efficient lookup of sums.
    fib_set = set(fib_sequence)

    # Step 2 & 3: Find and store unique valid combinations.
    # We use a set to store the sorted tuples of numbers in each valid combination
    # to automatically handle duplicates (e.g., a combination using the first '1'
    # vs. the same combination using the second '1').
    valid_combinations = set()
    min_len = 3
    max_len = 7

    for k in range(min_len, max_len + 1):
        # Generate combinations of size k from our list of Fibonacci numbers.
        for combo in itertools.combinations(fib_choices, k):
            # Calculate the sum of the current combination.
            s = sum(combo)

            # Check if the sum is a Fibonacci number.
            if s in fib_set:
                # If it is, add the sorted tuple of the combination to our set
                # to ensure uniqueness of the number set.
                valid_combinations.add(tuple(sorted(combo)))

    # Step 4: Print the results in a sorted and readable format.
    # Convert the set to a list and sort it for a consistent, ordered output.
    sorted_combinations = sorted(list(valid_combinations))

    for combo in sorted_combinations:
        # Format the numbers in the combination into an equation string.
        equation_str = " + ".join(map(str, combo))
        # Print the full equation.
        print(f"{equation_str} = {sum(combo)}")

    # Print the final count of unique combinations found.
    total_combinations = len(sorted_combinations)
    print(f"\nFound {total_combinations} combinations.")
    
    # Return the final answer in the specified format.
    print(f"<<<{total_combinations}>>>")

# Execute the function to find and print the combinations.
find_fibonacci_sum_combinations()