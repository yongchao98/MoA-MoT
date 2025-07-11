import itertools

def find_fibonacci_combinations():
    """
    This function finds and prints all combinations of 3 to 7 Fibonacci numbers
    whose sum is also a Fibonacci number, up to the 40th term.
    """
    fib_limit = 40
    
    # Step 1: Generate Fibonacci numbers.
    fibs = [1, 1]
    for i in range(2, fib_limit):
        fibs.append(fibs[i-1] + fibs[i-2])
    
    # Step 2 & 3: Create the pool for combinations and a set for fast lookup.
    # The pool includes F_1 to F_40. Note that fibs[0] is F_1 and fibs[39] is F_40.
    combination_pool = fibs[:fib_limit]
    fib_set = set(fibs)

    total_count = 0
    
    print(f"Searching for combinations of 3, 4, 5, 6, or 7 Fibonacci numbers (up to F_{fib_limit}) that sum to another Fibonacci number.")
    print("-" * 70)
    
    # Step 4: Iterate through combination sizes [3, 4, 5, 6, 7].
    for k in range(3, 8):
        # Generate all combinations of size k from the pool.
        # itertools.combinations handles the uniqueness of elements by their position,
        # so picking the two '1's is possible as they are at different positions.
        for combo in itertools.combinations(combination_pool, k):
            # Step 5: Calculate the sum and check if it's a Fibonacci number.
            current_sum = sum(combo)
            
            if current_sum in fib_set:
                total_count += 1
                # Format the numbers in the combo for printing.
                # The final equation should show each number.
                equation_parts = sorted(list(combo))
                equation_str = " + ".join(map(str, equation_parts))
                print(f"{equation_str} = {current_sum}")

    print("-" * 70)
    print(f"Total number of combinations found: {total_count}")

if __name__ == '__main__':
    find_fibonacci_combinations()