import itertools

def find_fibonacci_combinations():
    """
    Identifies and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number, up to the 40th Fibonacci number.
    It then prints the total count of such combinations.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The pool of numbers for combinations includes two '1's.
    n_limit = 40
    fib_pool = [1, 1]
    while len(fib_pool) < n_limit:
        next_fib = fib_pool[-1] + fib_pool[-2]
        fib_pool.append(next_fib)

    # Step 2: Create a set for efficient sum lookups.
    fib_set = set(fib_pool)
    max_fib_value = fib_pool[-1]

    total_count = 0
    found_equations = []

    print("Searching for combinations...")
    # Step 3: Iterate through combination sizes [3, 4, 5, 6, 7].
    for k in range(3, 8):
        # Step 4: Generate all combinations of size k from our pool.
        for combo in itertools.combinations(fib_pool, k):
            # Step 5: Calculate the sum and validate it.
            combo_sum = sum(combo)

            # Optimization: If the sum is already larger than the largest
            # possible Fibonacci number, we can skip it.
            if combo_sum > max_fib_value:
                continue

            # Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                total_count += 1
                # Format the equation string.
                equation = " + ".join(map(str, combo)) + f" = {combo_sum}"
                found_equations.append(equation)

    # Step 6: Print all the found valid equations.
    # Sorting makes the output deterministic and easier to read.
    if found_equations:
        print("\nFound the following valid combinations:\n")
        for eq in sorted(found_equations):
            print(eq)
    else:
        print("No valid combinations were found.")

    print(f"\nHow many such combinations are possible? \nTotal combinations: {total_count}")


if __name__ == "__main__":
    find_fibonacci_combinations()
<<<97>>>