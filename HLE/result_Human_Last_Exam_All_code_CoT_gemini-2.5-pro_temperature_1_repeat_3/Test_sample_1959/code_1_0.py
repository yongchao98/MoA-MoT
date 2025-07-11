import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of Fibonacci numbers (up to the 40th)
    of a specific length (3-7) whose sum is also a Fibonacci number.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with two 1s, which handles the duplicate '1' rule.
    n = 40
    fibs = [1, 1]
    while len(fibs) < n:
        fibs.append(fibs[-2] + fibs[-1])

    # Step 2: Create a set for efficient lookup of Fibonacci numbers.
    fib_set = set(fibs)

    # Initialize a counter for the valid combinations.
    found_combinations_count = 0
    
    # Define the constraints for the combination length.
    min_len = 3
    max_len = 7

    print("Searching for valid Fibonacci number combinations...")
    print("-" * 40)

    # Step 3: Iterate through the specified combination lengths.
    for k in range(min_len, max_len + 1):
        # Generate all combinations of length k from the Fibonacci list.
        # `itertools.combinations` handles the uniqueness of elements by their position.
        for combo in itertools.combinations(fibs, k):
            
            # Step 4: Calculate the sum and check if it's a Fibonacci number.
            combo_sum = sum(combo)
            
            if combo_sum in fib_set:
                # This is a valid combination. Increment the counter and print the result.
                found_combinations_count += 1
                
                # Format the numbers in the combination into a string like "1 + 1 + 2"
                equation_str = " + ".join(map(str, combo))
                
                print(f"{equation_str} = {combo_sum}")

    # Step 5: Print the final total count.
    print("-" * 40)
    print(f"Total number of combinations found: {found_combinations_count}")
    
    # Return the final answer in the required format.
    print(f"<<<{found_combinations_count}>>>")

solve_fibonacci_combinations()