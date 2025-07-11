import itertools

def find_fibonacci_combinations():
    """
    This function finds all combinations of 3 to 7 Fibonacci numbers
    (from the first 40) whose sum is also a Fibonacci number.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts F_1=1, F_2=1, ... up to F_40.
    n = 40
    fibs = [1, 1]
    while len(fibs) < n:
        fibs.append(fibs[-1] + fibs[-2])

    # For efficient lookup, create a set of the Fibonacci numbers.
    fib_set = set(fibs)

    # Initialize a counter for valid combinations.
    valid_combination_count = 0
    combination_lengths = [3, 4, 5, 6, 7]

    print(f"Finding combinations of Fibonacci numbers (up to the 40th term: {fibs[-1]})")
    print("The sum of these combinations must also be a Fibonacci number.\n")

    # Step 2: Iterate through the required combination sizes (3 to 7).
    for k in combination_lengths:
        # Step 3: Generate combinations.
        # itertools.combinations treats elements as unique based on their position.
        # This correctly handles the two '1's at the start of the sequence,
        # satisfying the condition that '1' can be used up to twice.
        for combo in itertools.combinations(fibs, k):
            
            # Step 4: Calculate the sum of the current combination.
            combo_sum = sum(combo)

            # Step 5: Check if the sum is a Fibonacci number.
            if combo_sum in fib_set:
                valid_combination_count += 1
                
                # Step 6: Print the valid combination and its sum.
                # Sorting the combo makes the output cleaner and consistent.
                equation_str = " + ".join(map(str, sorted(list(combo))))
                print(f"{equation_str} = {combo_sum}")

    # Final Step: Print the total count in the required format.
    print(f"\nFound a total of {valid_combination_count} combinations.")
    print(f"<<<{valid_combination_count}>>>")

# Execute the function to find and print the combinations and the final count.
find_fibonacci_combinations()