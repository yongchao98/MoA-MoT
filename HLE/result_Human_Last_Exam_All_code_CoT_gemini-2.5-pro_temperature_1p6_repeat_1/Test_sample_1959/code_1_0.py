import itertools

def find_fibonacci_sum_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number, using the first 40
    Fibonacci numbers as components.
    """

    # Step 1: Generate Fibonacci numbers for combinations (first 40)
    # The sequence starts with F_1=1, F_2=1.
    fib_components = []
    a, b = 1, 1
    for _ in range(40):
        fib_components.append(a)
        a, b = b, a + b

    # Step 2: Generate a larger set of Fibonacci numbers for checking the sum.
    # The maximum sum of a 7-number combo from the list is less than F_43.
    # We generate up to F_46 to be safe.
    fib_check_set = set()
    a, b = 1, 1
    for _ in range(46):
        fib_check_set.add(a)
        a, b = b, a + b

    count = 0
    print("Finding combinations of Fibonacci numbers that sum to a Fibonacci number...")
    print("-" * 70)

    # Step 3 & 4: Iterate through combination sizes, find, and check sums.
    for k in range(3, 8):  # Combinations of 3, 4, 5, 6, or 7 numbers
        for combo in itertools.combinations(fib_components, k):
            combo_sum = sum(combo)

            if combo_sum in fib_check_set:
                count += 1
                # Format and print the equation for the valid combination
                equation_str = " + ".join(map(str, combo))
                print(f"{equation_str} = {combo_sum}")

    # Step 5: Print the final count
    print("-" * 70)
    print(f"\nTotal number of combinations found: {count}")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()
<<<465>>>