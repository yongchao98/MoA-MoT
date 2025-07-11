import itertools

def find_fibonacci_sum_combinations():
    """
    Finds and prints all combinations of 3 to 7 unique Fibonacci numbers
    (with '1' allowed twice) that sum to another Fibonacci number, up to
    the 40th term in the sequence.
    """
    # Step 1: Generate Fibonacci numbers up to the 40th.
    # The list starts with two 1s to handle the duplication rule.
    fibs = [1, 1]
    a, b = 1, 1
    # Generate the sequence up to the 40th number.
    for _ in range(38):
        a, b = b, a + b
        fibs.append(b)

    # Step 2: Create a set for fast O(1) average time lookups.
    fibs_set = set(fibs)

    count = 0
    valid_combinations_output = []

    # Step 3 & 4: Iterate through combination lengths and generate combinations.
    for k in range(3, 8):
        # The `itertools.combinations` function will treat the two initial 1s as distinct
        # elements because they are at different positions in the input list.
        for combo in itertools.combinations(fibs, k):
            current_sum = sum(combo)

            # Step 5: Verify if the sum is a Fibonacci number.
            if current_sum in fibs_set:
                count += 1
                # Format the equation string for printing.
                equation = f"{' + '.join(map(str, combo))} = {current_sum}"
                valid_combinations_output.append(equation)

    # Step 6: Print all the valid equations.
    for line in valid_combinations_output:
        print(line)

    # Print the final count.
    print(f"\nFound {count} such combinations.")

if __name__ == '__main__':
    find_fibonacci_sum_combinations()