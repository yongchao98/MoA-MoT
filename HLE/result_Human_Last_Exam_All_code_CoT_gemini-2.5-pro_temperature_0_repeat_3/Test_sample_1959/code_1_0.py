import itertools

def find_fibonacci_sum_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number, up to the 40th Fibonacci number.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with F_1=1, F_2=1. This list naturally allows '1'
    # to be used up to twice in combinations.
    fib_numbers = [1, 1]
    while len(fib_numbers) < 40:
        next_fib = fib_numbers[-1] + fib_numbers[-2]
        fib_numbers.append(next_fib)

    # Create a set for fast lookups to check if a sum is a Fibonacci number.
    fib_set = set(fib_numbers)

    count = 0
    print("Finding combinations of Fibonacci numbers whose sum is also a Fibonacci number...")
    print("-" * 70)

    # Step 2: Iterate through combination lengths from 3 to 7.
    for k in range(3, 8):
        # Generate all unique combinations of length k.
        for combo in itertools.combinations(fib_numbers, k):
            # Step 3: Calculate the sum of the current combination.
            current_sum = sum(combo)

            # Step 4: Check if the sum is a Fibonacci number.
            if current_sum in fib_set:
                # If it is, increment the count and print the result.
                count += 1
                # Format the output string as an equation.
                # The problem asks to output each number in the final equation.
                equation_str = " + ".join(map(str, combo)) + f" = {current_sum}"
                print(equation_str)

    # Step 5: Print the total count.
    print("-" * 70)
    print(f"\nTotal number of combinations found: {count}")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()
<<<48>>>