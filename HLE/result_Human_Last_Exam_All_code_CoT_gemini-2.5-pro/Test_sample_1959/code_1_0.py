import itertools

# Plan:
# 1. Generate Fibonacci numbers: Create a list of the first 40 Fibonacci numbers, starting with 1, 1. This list will be the source for our combinations.
# 2. Create a Fibonacci set: Convert the list of Fibonacci numbers into a set for fast O(1) average time complexity lookups. This will be used to check if the sum of a combination is a Fibonacci number.
# 3. Iterate through combination sizes: Loop through the required number of elements in a combination, which are 3, 4, 5, 6, and 7.
# 4. Generate combinations: For each size, use `itertools.combinations` to generate all possible unique combinations from our list of 40 Fibonacci numbers. `itertools.combinations` handles the non-duplication rule automatically. Because our source list contains two '1's, combinations can include up to two '1's.
# 5. Check the sum: For each generated combination, calculate its sum.
# 6. Validate the sum: Check if the calculated sum exists in the Fibonacci set created in step 2.
# 7. Print results: If the sum is a valid Fibonacci number, print the combination as a formatted equation (e.g., "1 + 2 + 5 = 8").
# 8. Count and report: Keep a running total of all valid combinations found. After checking all possibilities, print the final count.

# Step 1: Generate the first 40 Fibonacci numbers
n = 40
fib_numbers = [1, 1]
# We start the loop from 2 because the first two numbers (at index 0 and 1) are already defined.
for i in range(2, n):
    next_fib = fib_numbers[i-1] + fib_numbers[i-2]
    fib_numbers.append(next_fib)

# Step 2: Create a set for efficient sum checking
fib_set = set(fib_numbers)

found_combinations_count = 0
print("Finding all combinations of 3 to 7 Fibonacci numbers whose sum is also a Fibonacci number...")
print("-" * 80)

# Step 3: Iterate through combination sizes (k=3 to 7)
for k in range(3, 8):
    # Step 4: Generate combinations of size k from the list of Fibonacci numbers
    # The `itertools.combinations` function works on the positions of elements in the list.
    # Since `fib_numbers` starts with [1, 1, ...], this correctly allows '1' to be used twice.
    for combo in itertools.combinations(fib_numbers, k):
        # Step 5: Calculate the sum of the current combination
        current_sum = sum(combo)

        # Step 6: Validate the sum against the set of Fibonacci numbers
        if current_sum in fib_set:
            # Step 7: Increment the counter and print the formatted equation
            found_combinations_count += 1
            # Sort the numbers in the combination for a clean, consistent output format
            sorted_combo = sorted(list(combo))
            # Output each number in the final equation as requested
            equation_str = " + ".join(map(str, sorted_combo)) + f" = {current_sum}"
            print(equation_str)

# Step 8: Print the final total count
print("-" * 80)
print(f"Total number of combinations found: {found_combinations_count}")