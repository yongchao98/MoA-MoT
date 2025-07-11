import itertools

def find_fibonacci_sum_combinations():
    """
    Identifies and prints all combinations of 3 to 7 Fibonacci numbers
    (from the first 40) whose sum is also a Fibonacci number. It then
    prints the total count of such combinations.
    """
    # Step 1: Generate Fibonacci numbers.
    # We generate up to F_46 to ensure our check-set is large enough for the sums.
    fib_sequence = [1, 1]
    while len(fib_sequence) < 46:
        next_fib = fib_sequence[-1] + fib_sequence[-2]
        fib_sequence.append(next_fib)

    # The pool of numbers for our combinations is the first 40 Fibonacci numbers.
    # This list starts with [1, 1, 2, 3, 5, ...] up to F_40.
    combination_pool = fib_sequence[:40]

    # A set for efficient checking of whether a number is in the Fibonacci sequence.
    fib_check_set = set(fib_sequence)

    found_combinations_count = 0

    print("Found the following combinations:")

    # Step 2 & 3: Iterate through combination lengths, generate, and validate.
    for k in range(3, 8):
        for combo in itertools.combinations(combination_pool, k):
            current_sum = sum(combo)

            # Check if the sum is a Fibonacci number.
            if current_sum in fib_check_set:
                found_combinations_count += 1
                # As required, print each number in the final equation.
                equation_str = " + ".join(map(str, combo))
                print(f"{equation_str} = {current_sum}")

    # Step 4: Output the final count.
    print(f"\nTotal number of combinations found: {found_combinations_count}")


if __name__ == '__main__':
    find_fibonacci_sum_combinations()