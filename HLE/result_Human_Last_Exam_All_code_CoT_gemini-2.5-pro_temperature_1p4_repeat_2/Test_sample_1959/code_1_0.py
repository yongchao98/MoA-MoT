import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers
    (from the first 40) whose sum is also a Fibonacci number.
    """
    # Step 1: Generate Fibonacci numbers.
    # The max sum of the 7 largest numbers from the first 40 is less than F_42.
    # So generating up to F_42 is sufficient for the lookup set.
    def generate_fibs(n):
        fibs = [1, 1]
        a, b = 1, 1
        for _ in range(n - 2):
            a, b = b, a + b
            fibs.append(b)
        return fibs

    all_fibs = generate_fibs(42)
    fib_set = set(all_fibs)

    # Step 2: Create the pool of Fibonacci numbers for combinations.
    # The first 40 Fibonacci numbers (F_1 to F_40). This list includes '1' twice.
    fib_source = all_fibs[:40]

    # Step 3: Find combinations, check sum, and count.
    count = 0
    found_combinations = []

    # Loop through combination sizes from 3 to 7.
    for k in range(3, 8):
        # Generate all combinations of size k from the source pool.
        for combo in itertools.combinations(fib_source, k):
            # Calculate the sum of the combination.
            combo_sum = sum(combo)

            # Step 4: Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                count += 1
                # Sort the combination for consistent output formatting.
                sorted_combo = sorted(list(combo))
                found_combinations.append((sorted_combo, combo_sum))

    # Step 5: Print the results.
    # Sort the final list of combinations for a deterministic output.
    # Sort first by sum, then by length of combo, then by the combo itself.
    found_combinations.sort(key=lambda x: (x[1], len(x[0]), x[0]))

    for combo, s in found_combinations:
        # Create the equation string and print it.
        # Each number in the final equation is printed as requested.
        equation_str = " + ".join(map(str, combo)) + f" = {s}"
        print(equation_str)

    print(f"\nTotal number of combinations found: {count}")
    
    # Return the final count for the specified answer format
    return count

if __name__ == '__main__':
    # The final answer is the total count, which we will wrap in <<<>>>
    final_count = solve_fibonacci_combinations()
    # The following line is commented out to prevent printing the answer tag
    # directly into the console output when the script is run.
    # print(f"<<<{final_count}>>>")
