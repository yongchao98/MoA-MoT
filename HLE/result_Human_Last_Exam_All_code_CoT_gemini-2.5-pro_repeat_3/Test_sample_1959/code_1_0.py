import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers whose sum
    is also a Fibonacci number, up to the 40th Fibonacci number.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1, which handles the constraint that
    # '1' can be used twice in a combination.
    n = 40
    fib_numbers = [1, 1]
    for i in range(2, n):
        fib_numbers.append(fib_numbers[i-1] + fib_numbers[i-2])

    # Step 2: Create a set for efficient sum verification.
    fib_set = set(fib_numbers)

    # Initialize a list to store the results and a counter.
    found_equations = []
    total_count = 0

    # Step 3: Iterate through combination lengths from 3 to 7.
    for k in range(3, 8):
        # Step 4: Generate all unique combinations of size k.
        # itertools.combinations handles the uniqueness constraint automatically.
        for combo in itertools.combinations(fib_numbers, k):
            combo_sum = sum(combo)

            # Step 5: Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                total_count += 1
                # Format the equation string with numbers sorted for consistency.
                equation_str = " + ".join(map(str, sorted(list(combo)))) + f" = {combo_sum}"
                found_equations.append(equation_str)

    # Step 6: Print all the discovered equations.
    if found_equations:
        print("Found the following combinations:")
        for eq in found_equations:
            print(eq)
    else:
        print("No combinations were found matching the criteria.")

    # Print the final count.
    print(f"\nHow many such combinations are possible? {total_count}")

if __name__ == "__main__":
    solve_fibonacci_combinations()