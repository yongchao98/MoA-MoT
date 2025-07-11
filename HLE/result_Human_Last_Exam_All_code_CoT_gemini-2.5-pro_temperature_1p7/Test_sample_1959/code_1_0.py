import itertools

def solve_fibonacci_combinations():
    """
    This function finds and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number, up to the 40th Fibonacci number.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1, and the max number allowed is F_40.
    n = 40
    fib_numbers = [1, 1]
    # We stop when the list has 40 numbers.
    while len(fib_numbers) < n:
        next_fib = fib_numbers[-1] + fib_numbers[-2]
        fib_numbers.append(next_fib)

    # Step 2: Create a set for efficient O(1) average time lookups.
    fib_set = set(fib_numbers)

    # Initialize a counter for valid combinations.
    valid_combination_count = 0

    print("Searching for combinations... Found the following:\n")

    # Step 3: Iterate through combination sizes (3, 4, 5, 6, 7).
    for k in range(3, 8):
        # Step 4: Generate combinations of size k.
        # The list 'fib_numbers' correctly handles the duplicate '1'.
        for combo in itertools.combinations(fib_numbers, k):
            # Step 5: Calculate the sum of the current combination.
            combo_sum = sum(combo)
            
            # Step 6: Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                valid_combination_count += 1
                # Format the numbers in the combination for printing.
                equation_str = " + ".join(map(str, combo))
                # Step 7: Print the resulting equation.
                print(f"{equation_str} = {combo_sum}")

    # After checking all combinations, print the final count.
    print(f"\nTotal number of combinations found: {valid_combination_count}")

if __name__ == '__main__':
    solve_fibonacci_combinations()