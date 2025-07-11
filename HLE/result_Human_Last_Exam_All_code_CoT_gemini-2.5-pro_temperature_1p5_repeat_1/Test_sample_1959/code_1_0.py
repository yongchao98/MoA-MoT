import itertools

def find_fibonacci_sum_combinations():
    """
    This script finds all combinations of 3 to 7 distinct Fibonacci numbers
    (where 1 can be used twice) whose sum is also a Fibonacci number.
    The numbers in the combination are no larger than the 40th Fibonacci number.
    """
    
    # Step 1: Generate the pool of Fibonacci numbers to choose from (F_1 to F_40).
    # The sequence is 1, 1, 2, 3, 5, ..., F(40).
    fib_pool = [1, 1]
    a, b = 1, 1
    # Generate F_3 through F_40
    for _ in range(38):
        c = a + b
        fib_pool.append(c)
        a, b = b, c

    # Step 2: Generate a set of Fibonacci numbers for quick sum lookups.
    # The maximum possible sum is sum(fib_pool[-7:]), so we need a set that covers that.
    fib_check_set = set()
    x, y = 1, 1
    # Calculate an upper bound for the sum to ensure the check set is large enough.
    max_sum_possible = sum(sorted(fib_pool, reverse=True)[:7])
    while x <= max_sum_possible:
        fib_check_set.add(x)
        x, y = y, x + y
    # The loop might exit before the next fibonacci number (the sum itself) is added.
    fib_check_set.add(x)


    # Step 3 & 5: Store unique, valid combinations in a set.
    # We use a set of sorted tuples to ensure uniqueness, as combinations of values
    # like {1, 2, 3} might be generated multiple times due to the two '1's.
    found_combinations = set()

    # Step 4: Iterate through combination sizes from 3 to 7.
    for k in range(3, 8):
        # Generate all combinations of size k from the pool.
        for combo in itertools.combinations(fib_pool, k):
            
            # Calculate the sum of the current combination.
            current_sum = sum(combo)

            # Check if the sum is a Fibonacci number.
            if current_sum in fib_check_set:
                # If it is, add the sorted tuple representation of the combination
                # to our set of results to ensure uniqueness.
                sorted_combo = tuple(sorted(combo))
                found_combinations.add(sorted_combo)

    # Step 6: Print the results in a sorted order.
    # First, convert the set to a list to be able to sort it.
    sorted_results = sorted(list(found_combinations))

    print("Found the following combinations:\n")
    for combo in sorted_results:
        combo_sum = sum(combo)
        # Remember to output each number in the final equation
        equation = " + ".join(map(str, combo)) + f" = {combo_sum}"
        print(equation)

    # Print the final count.
    total_count = len(found_combinations)
    print(f"\nHow many such combinations are possible? {total_count}")


if __name__ == '__main__':
    find_fibonacci_sum_combinations()
<<<51>>>