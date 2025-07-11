import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers (from the first 40)
    that sum to another Fibonacci number.
    """
    # Step 1: Generate the pool of Fibonacci numbers up to the 40th.
    # The problem starts the sequence with F(1)=1, F(2)=1.
    # The list will naturally contain two 1s.
    fibs_pool = [1, 1]
    while len(fibs_pool) < 40:
        next_fib = fibs_pool[-1] + fibs_pool[-2]
        fibs_pool.append(next_fib)

    # Step 2: Create a larger set of Fibonacci numbers for efficient sum checking.
    # The maximum sum could be larger than the 40th Fibonacci number. We extend
    # the sequence to a safe upper bound to ensure all valid sums are included.
    # The sum of the 7 largest numbers in the pool is a safe upper bound.
    max_possible_sum = sum(sorted(fibs_pool, reverse=True)[:7])
    
    fib_check_list = list(fibs_pool)
    while fib_check_list[-1] <= max_possible_sum:
        fib_check_list.append(fib_check_list[-1] + fib_check_list[-2])
    fib_check_set = set(fib_check_list)

    total_count = 0
    
    print("Finding all combinations of 3 to 7 Fibonacci numbers whose sum is also a Fibonacci number...")
    print("The pool of numbers is the first 40 Fibonacci numbers (up to 102,334,155).")
    print("-" * 70)

    # Step 3 & 4: Iterate through combination lengths, generate combinations, and validate them.
    for k in range(3, 8):
        # itertools.combinations treats elements as unique based on their position,
        # so it correctly handles the two '1's in the fibs_pool.
        for combo in itertools.combinations(fibs_pool, k):
            combo_sum = sum(combo)
            
            # Check if the sum is a Fibonacci number.
            if combo_sum in fib_check_set:
                total_count += 1
                # Format the numbers in the combo for printing the equation.
                equation_str = " + ".join(map(str, sorted(list(combo))))
                print(f"{equation_str} = {combo_sum}")

    # Step 5: Print the final total count.
    print("-" * 70)
    print(f"Total number of combinations found: {total_count}")
    
    # Return the final answer in the specified format.
    print(f"\n<<<{total_count}>>>")


solve_fibonacci_combinations()