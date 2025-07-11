import itertools

def find_fibonacci_sum_combinations():
    """
    Identifies all combinations of 3 to 7 Fibonacci numbers whose sum is also a Fibonacci number.
    The first 40 Fibonacci numbers are used as the pool for the combinations.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1 as per the problem description.
    fibs = [1, 1]
    while len(fibs) < 40:
        next_fib = fibs[-1] + fibs[-2]
        fibs.append(next_fib)

    # Step 2: Create a set for fast lookup of Fibonacci numbers.
    fib_set = set(fibs)

    # This set will store the unique combinations that meet the criteria.
    # Using a set of tuples prevents duplicates, e.g., (1, 2, 5) and (2, 1, 5).
    valid_combinations_set = set()

    # Step 3 & 4: Iterate through combination sizes (3 to 7) and generate combinations.
    for k in range(3, 8):
        for combo in itertools.combinations(fibs, k):
            # Step 5: Check if the sum of the combination is a Fibonacci number.
            current_sum = sum(combo)
            if current_sum in fib_set:
                # Add the combination as a sorted tuple to the set to ensure uniqueness.
                valid_combinations_set.add(combo)

    # Sort the results for a consistent and readable output.
    sorted_combinations = sorted(list(valid_combinations_set))

    # Step 6: Print each valid combination in the required format.
    for combo in sorted_combinations:
        # Each number in the combo needs to be part of the final equation string.
        equation_parts = [str(n) for n in combo]
        equation = " + ".join(equation_parts)
        print(f"{equation} = {sum(combo)}")

    # Finally, print the total count of unique combinations found.
    total_count = len(sorted_combinations)
    print(f"\nTotal number of combinations found: {total_count}")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()