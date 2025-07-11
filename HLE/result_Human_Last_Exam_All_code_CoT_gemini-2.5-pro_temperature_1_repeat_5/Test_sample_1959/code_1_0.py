import itertools

def find_fibonacci_sum_combinations():
    """
    This script finds all combinations of 3 to 7 Fibonacci numbers whose sum
    is also a Fibonacci number, based on the problem's constraints.
    """
    # 1. Generate the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1, 2, 3, ...
    n = 40
    fib_numbers = [1, 1]
    while len(fib_numbers) < n:
        next_fib = fib_numbers[-1] + fib_numbers[-2]
        fib_numbers.append(next_fib)

    # 2. Create a set of Fibonacci numbers for efficient lookups.
    fib_set = set(fib_numbers)

    total_count = 0
    
    print(f"Finding all combinations of 3, 4, 5, 6, or 7 Fibonacci numbers (from the first {n}) that sum to a Fibonacci number.")
    print("-" * 80)

    # 3. Iterate through combination lengths (3 to 7).
    for k in range(3, 8):
        # 4. Generate all combinations of length k from the list of Fibonacci numbers.
        # The list [1, 1, 2, ...] ensures '1' can be used twice, while other numbers
        # are used at most once per combination.
        combinations = itertools.combinations(fib_numbers, k)
        
        for combo in combinations:
            # 5. Calculate the sum of the current combination.
            combo_sum = sum(combo)

            # 6. Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                total_count += 1
                # Format and print the equation for the valid combination.
                # Sorting the combo makes the output easier to read.
                equation_str = " + ".join(map(str, sorted(list(combo)))) + f" = {combo_sum}"
                print(equation_str)

    # 7. Print the final total count.
    print("-" * 80)
    print(f"Total number of such combinations found: {total_count}")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()