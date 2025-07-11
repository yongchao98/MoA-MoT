import itertools

def solve_fibonacci_combinations():
    """
    Finds and counts all combinations of 3 to 7 Fibonacci numbers whose sum is
    also a Fibonacci number, up to the 40th term.
    """
    # Step 1: Generate the first 40 Fibonacci numbers.
    # The sequence starts with 1, 1, so we handle the rule that '1' can be used twice.
    fibs = [1, 1]
    while len(fibs) < 40:
        fibs.append(fibs[-1] + fibs[-2])

    # Step 2: Create a set for efficient O(1) average time lookups.
    fib_set = set(fibs)

    count = 0
    
    # Step 3: Iterate through combination lengths from 3 to 7.
    for k in range(3, 8):
        # Generate all combinations of size k from the list of Fibonacci numbers.
        for combo in itertools.combinations(fibs, k):
            # Step 4: Calculate the sum and check if it's a Fibonacci number.
            combo_sum = sum(combo)
            if combo_sum in fib_set:
                # Step 5: If it is, increment the count and print the equation.
                count += 1
                # The final equation requires printing each number.
                equation_parts = [str(num) for num in combo]
                print(f"{' + '.join(equation_parts)} = {combo_sum}")

    # Step 6: Print the total count.
    print(f"\nTotal number of combinations found: {count}")
    
    # The final answer format as requested by the user.
    print(f"\n<<<{count}>>>")

if __name__ == '__main__':
    solve_fibonacci_combinations()