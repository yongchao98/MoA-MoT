import itertools

def find_fibonacci_sum_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number.
    The Fibonacci numbers are capped at the 40th term.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with F_1=1, F_2=1.
    fibs = [1, 1]
    while len(fibs) < 40:
        fibs.append(fibs[-1] + fibs[-2])

    # Create a set for fast lookups to check if a number is in the sequence.
    fib_set = set(fibs)

    count = 0
    print("Found the following combinations:")

    # Step 2: Iterate through combination lengths (3, 4, 5, 6, 7).
    for k in range(3, 8):
        # Step 3: Generate combinations of length k from the Fibonacci numbers.
        # itertools.combinations handles the duplicate '1's correctly.
        for combo in itertools.combinations(fibs, k):
            # Step 4: Calculate the sum of the current combination.
            current_sum = sum(combo)

            # Step 5: Check if the sum is a Fibonacci number.
            if current_sum in fib_set:
                # Step 6: If it is, print the equation and increment the counter.
                count += 1
                # Create the equation string.
                equation_str = " + ".join(map(str, combo)) + f" = {current_sum}"
                print(equation_str)

    # Step 7: Print the total count of valid combinations found.
    print(f"\nTotal number of combinations found: {count}")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()