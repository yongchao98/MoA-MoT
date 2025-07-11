import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers whose sum is also a Fibonacci number.
    """
    # Step 1: Generate the source list of the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1, as per the problem description.
    # This list will be used to generate the combinations.
    fib_source_list = [1, 1]
    while len(fib_source_list) < 40:
        next_fib = fib_source_list[-1] + fib_source_list[-2]
        fib_source_list.append(next_fib)

    # Step 2: Generate a larger set of Fibonacci numbers for efficient sum checking.
    # The largest possible sum will be the sum of the 7 largest numbers from our source list.
    # We generate Fibonacci numbers up to a value safely greater than that potential sum.
    fib_check_list = list(fib_source_list)
    # The sum of the last 7 elements is < 7*F_40 ~ 7.1e8. F_45 is > 1.1e9, which is a safe upper bound.
    while fib_check_list[-1] < 1_200_000_000:
        next_fib = fib_check_list[-1] + fib_check_list[-2]
        fib_check_list.append(next_fib)
    fib_check_set = set(fib_check_list)

    # Step 3 & 4: Iterate through combination lengths, generate combinations, and check their sum.
    found_combinations_count = 0
    solution_equations = []
    combination_lengths = range(3, 8) # Combination lengths of 3, 4, 5, 6, or 7

    for k in combination_lengths:
        for combo in itertools.combinations(fib_source_list, k):
            combo_sum = sum(combo)

            # Check if the sum is a Fibonacci number
            if combo_sum in fib_check_set:
                found_combinations_count += 1
                
                # Sort numbers within the combination for consistent output
                sorted_combo = sorted(list(combo))
                
                # Format the equation string
                equation = " + ".join(map(str, sorted_combo)) + f" = {combo_sum}"
                solution_equations.append(equation)
    
    # Step 5: Print the results.
    # Sorting the final list of equations makes the output tidy and easy to inspect.
    solution_equations.sort()
    for eq in solution_equations:
        print(eq)

    print(f"\nTotal number of combinations found: {found_combinations_count}")
    
    # Return the final count in the specified format
    print(f"<<<{found_combinations_count}>>>")

if __name__ == "__main__":
    solve_fibonacci_combinations()
