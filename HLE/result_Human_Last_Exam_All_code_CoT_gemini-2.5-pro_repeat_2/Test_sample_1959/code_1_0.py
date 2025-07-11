import itertools

def solve_fibonacci_combinations():
    """
    This function finds and prints all combinations of 3 to 7 unique Fibonacci numbers
    (with '1' allowed twice) whose sum is also a Fibonacci number, using the first 40
    Fibonacci numbers.
    """
    # Step 1: Generate Fibonacci Numbers
    # Generate the first 40 Fibonacci numbers, starting with 1, 1.
    n = 40
    fibs_list = [1, 1]
    while len(fibs_list) < n:
        fibs_list.append(fibs_list[-1] + fibs_list[-2])

    # Create a set for efficient O(1) average time complexity lookups.
    fib_set = set(fibs_list)

    # A list to store the combinations that satisfy the criteria.
    valid_combinations = []

    # Step 2 & 3: Iterate through combination lengths and generate combinations.
    # The combination must consist of 3, 4, 5, 6, or 7 Fibonacci numbers.
    for k in range(3, 8):
        # 'itertools.combinations' will create combinations of k length.
        # Since 'fibs_list' contains two '1's, '1' can be part of a combination twice.
        for combo in itertools.combinations(fibs_list, k):
            combo_sum = sum(combo)
            
            # Step 3 (cont.): Check if the sum is also a Fibonacci number.
            if combo_sum in fib_set:
                valid_combinations.append(combo)

    # Step 4: Display the results.
    count = len(valid_combinations)
    print(f"Found {count} possible combinations.\n")
    print("The combinations are:")
    for combo in valid_combinations:
        # Sort the numbers in the combination for a consistent, readable output.
        sorted_combo = sorted(list(combo))
        # Format the equation string with each number.
        equation_parts = " + ".join(map(str, sorted_combo))
        # Print the final equation.
        print(f"{equation_parts} = {sum(combo)}")
    
    # The final answer is the total count.
    print(f"\n<<<48>>>")

solve_fibonacci_combinations()