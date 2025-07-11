import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers (up to F(40))
    that sum to another Fibonacci number. Finally, prints the total count.
    """

    def generate_fibs(n):
        """
        Generates the first n Fibonacci numbers, using the definition F(1)=1, F(2)=1.
        """
        if n <= 0:
            return []
        fibs = [1, 1]
        while len(fibs) < n:
            fibs.append(fibs[-1] + fibs[-2])
        return fibs[:n]

    # 1. Generate the pool of Fibonacci numbers for making combinations.
    #    The pool consists of the first 40 Fibonacci numbers.
    fib_pool = generate_fibs(40)

    # 2. Generate a larger set of Fibonacci numbers for checking sums.
    #    The max sum of 7 numbers from the pool is < F(42). F(45) is a safe upper bound.
    fib_check_set = set(generate_fibs(45))

    # 3. Initialize a counter for valid combinations.
    valid_combinations_count = 0

    print("Finding combinations of Fibonacci numbers that sum to a Fibonacci number...")

    # 4. Iterate through combination sizes (k=3 to 7).
    for k in range(3, 8):
        # 5. Generate all unique combinations of size k from the pool.
        #    itertools.combinations treats the two '1's as distinct, allowing
        #    them to be chosen together in a combination.
        for combo in itertools.combinations(fib_pool, k):
            current_sum = sum(combo)

            # 6. Check if the sum is a Fibonacci number (and not in the combo itself to avoid trivial sums like 3+5=8 if 8 is in combo).
            #    However, the problem does not forbid the sum from being one of the numbers.
            #    A simple check `current_sum in fib_check_set` is sufficient.
            if current_sum in fib_check_set:
                valid_combinations_count += 1
                # 7. Print the found combination and its sum.
                equation_str = " + ".join(map(str, combo))
                print(f"{equation_str} = {current_sum}")

    # 8. Print the final count of valid combinations.
    print(f"\nTotal number of combinations found: {valid_combinations_count}")
    print(f"<<<{valid_combinations_count}>>>")

# Execute the solution
solve_fibonacci_combinations()