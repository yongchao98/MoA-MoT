import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers whose sum
    is also a Fibonacci number.
    """
    # Step 1: Generate Fibonacci numbers.
    # We need the first 40 for the pool. We need more for the check set.
    # The max sum of 7 largest numbers from the F(40) pool is ~258M (F(42)).
    # We'll generate up to F(45) to be safe.
    fib_list = [1, 1]
    while len(fib_list) < 45:
        next_fib = fib_list[-1] + fib_list[-2]
        fib_list.append(next_fib)

    # The pool for combinations is the first 40 Fibonacci numbers.
    fib_pool = fib_list[:40]

    # The check set is for quick lookups to see if a sum is a Fibonacci number.
    fib_check_set = set(fib_list)

    # Use a set to store found combinations to automatically handle duplicates
    # that arise from the two '1's in the pool.
    found_combinations_set = set()

    # Step 2 & 3: Iterate through combination sizes, generate, and check.
    min_combo_size = 3
    max_combo_size = 7

    for k in range(min_combo_size, max_combo_size + 1):
        for combo in itertools.combinations(fib_pool, k):
            combo_sum = sum(combo)

            if combo_sum in fib_check_set:
                # Store the sorted tuple to ensure uniqueness of the values.
                sorted_combo = tuple(sorted(combo))
                found_combinations_set.add(sorted_combo)

    # Step 4: Print the results.
    # Sort the final list of tuples for a consistent, ordered output.
    sorted_results = sorted(list(found_combinations_set), key=lambda x: (len(x), sum(x), x))
    
    for combo in sorted_results:
        # Construct and print the equation for each valid combination.
        # This fulfills the requirement to "output each number in the final equation".
        equation_str = " + ".join(map(str, combo))
        combo_sum = sum(combo)
        print(f"{equation_str} = {combo_sum}")

    # Print the final count.
    final_count = len(found_combinations_set)
    print(f"\nFound {final_count} unique combinations.")

if __name__ == '__main__':
    solve_fibonacci_combinations()
