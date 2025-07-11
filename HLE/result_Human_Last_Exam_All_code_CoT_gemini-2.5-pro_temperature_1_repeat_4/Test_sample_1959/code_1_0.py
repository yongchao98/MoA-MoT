import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 unique Fibonacci numbers
    (from the first 40) whose sum is also a Fibonacci number.
    """
    # Step 1: Generate Fibonacci numbers
    # The pool for making combinations will be the first 40 Fibonacci numbers.
    max_n_in_combo = 40
    # The set for checking the sum needs to be larger. The sum of the 7 largest
    # numbers from the pool of 40 is less than the 47th Fibonacci number.
    max_n_for_check = 47

    fib_pool = []
    fib_list_for_check = []
    a, b = 1, 1
    for i in range(max_n_for_check):
        # The first 40 numbers (including two 1s) go into our combination pool
        if i < max_n_in_combo:
            fib_pool.append(a)
        # A larger list is created for the check set
        fib_list_for_check.append(a)
        a, b = b, a + b

    # Using a set for the check provides fast O(1) average time complexity for lookups
    fib_check_set = set(fib_list_for_check)

    # Step 2 & 3: Iterate through combination sizes, find, and print solutions
    found_combinations_count = 0
    min_combo_size = 3
    max_combo_size = 7

    print("--- Found Combinations ---")

    for k in range(min_combo_size, max_combo_size + 1):
        # Generate all unique combinations of size k from the pool
        for combo in itertools.combinations(fib_pool, k):
            current_sum = sum(combo)

            # Check if the sum is a Fibonacci number
            if current_sum in fib_check_set:
                found_combinations_count += 1
                # Format the numbers in the combination into a string like "1 + 1 + 2"
                equation_parts = " + ".join(map(str, combo))
                print(f"{equation_parts} = {current_sum}")

    # Step 4: Final Output
    print("\n--- Summary ---")
    print(f"The maximum Fibonacci number allowed in a combination was: {fib_pool[-1]:,}")
    print(f"Combinations were created using the first {len(fib_pool)} Fibonacci numbers.")
    print(f"Found a total of {found_combinations_count} combinations.")

    # Final answer in the required format
    print(f"<<<{found_combinations_count}>>>")

if __name__ == '__main__':
    solve_fibonacci_combinations()