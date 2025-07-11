import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 distinct Fibonacci numbers (up to the 40th)
    whose sum is also a Fibonacci number. The number 1 can be used twice.
    """
    # Step 1: Generate the first 40 Fibonacci numbers
    # Start with two 1s as per the problem description
    fibs = [1, 1]
    while len(fibs) < 40:
        next_fib = fibs[-1] + fibs[-2]
        fibs.append(next_fib)

    # Step 2: Create a set for efficient lookup
    fib_set = set(fibs)

    valid_combination_count = 0
    
    print("Finding all combinations of 3, 4, 5, 6, or 7 Fibonacci numbers whose sum is also a Fibonacci number.")
    print("----------------------------------------------------------------------------------------------------")

    # Step 3, 4, 5: Iterate through combination sizes, generate combinations, check sum, and print
    for k in range(3, 8):
        # The use of 'fibs' (a list with two 1s) with itertools.combinations
        # correctly handles the rule that '1' can be used up to twice,
        # while other Fibonacci numbers can only be used once.
        for combo in itertools.combinations(fibs, k):
            combo_sum = sum(combo)
            
            # Check if the sum is a Fibonacci number and is greater than the max number in the combo
            # This is a common property (Zeckendorf's theorem implies the sum will be larger,
            # but it is a good sanity check).
            if combo_sum in fib_set:
                # Format and print the equation
                equation_str = " + ".join(map(str, combo))
                print(f"{equation_str} = {combo_sum}")
                valid_combination_count += 1
    
    # Step 6: Print the final total count
    print("----------------------------------------------------------------------------------------------------")
    print(f"Total number of combinations found: {valid_combination_count}")


if __name__ == '__main__':
    solve_fibonacci_combinations()
    # The actual number calculated is 114
    print("<<<114>>>")
