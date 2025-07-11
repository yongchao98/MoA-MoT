import itertools

def find_fibonacci_sum_combinations():
    """
    This script identifies and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number.

    The search space is limited to the first 40 Fibonacci numbers.
    """
    print("Initializing...")

    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts [1, 1, ...], which allows the number 1 to be used twice
    # as per the problem's constraints, since itertools.combinations treats elements
    # as unique based on their position in the input list.
    n = 40
    fib_pool = [1, 1]
    a, b = 1, 1
    # We already have the first two, so generate the next 38.
    for _ in range(n - 2):
        a, b = b, a + b
        fib_pool.append(b)

    # Step 2: Create a set for fast lookup of Fibonacci numbers.
    fib_set = set(fib_pool)
    
    combination_count = 0
    print("Searching for valid combinations. This may take a moment...")
    print("-" * 30)

    # Step 4: Iterate through combination lengths from 3 to 7.
    for k in range(3, 8):
        # Generate all combinations of length k from our pool of Fibonacci numbers.
        for combo in itertools.combinations(fib_pool, k):
            current_sum = sum(combo)

            # Step 5: Check if the sum is a Fibonacci number.
            if current_sum in fib_set:
                combination_count += 1
                # Format the numbers in the combination for printing.
                combo_str = " + ".join(map(str, sorted(list(combo))))
                print(f"{combo_str} = {current_sum}")

    # Step 6: Print the final total count.
    print("-" * 30)
    print(f"\nTotal number of combinations found: {combination_count}")


if __name__ == "__main__":
    find_fibonacci_sum_combinations()
