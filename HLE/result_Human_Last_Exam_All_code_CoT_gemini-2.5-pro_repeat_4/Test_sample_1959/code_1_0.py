import itertools

def find_fibonacci_sum_combinations():
    """
    This function finds and prints all combinations of 3 to 7 Fibonacci numbers
    (from the first 40) whose sum is also a Fibonacci number.
    It then prints the total count of such combinations.
    """
    # Step 1: Generate Fibonacci numbers.
    # We need the first 40 for creating combinations.
    # The sum can be larger, so we generate up to F_45 for the lookup set.
    fibs = [1, 1]
    while len(fibs) < 45:
        fibs.append(fibs[-1] + fibs[-2])

    # The first 40 Fibonacci numbers are allowed in the combination.
    combination_fibs = fibs[:40]
    # Use a set for efficient lookup of sums.
    fib_lookup_set = set(fibs)

    found_combinations = []
    
    # Step 2 & 3: Iterate through combination sizes and generate combinations.
    for k in range(3, 8):
        for combo in itertools.combinations(combination_fibs, k):
            # Step 4: Check if the sum is a Fibonacci number.
            combo_sum = sum(combo)
            if combo_sum in fib_lookup_set:
                found_combinations.append((combo, combo_sum))

    # Step 5: Print the found combinations and their sums.
    print("Found the following combinations:")
    for combo, combo_sum in found_combinations:
        # Format the equation string.
        equation_str = " + ".join(map(str, combo))
        print(f"{equation_str} = {combo_sum}")

    # Step 6: Print the total count.
    count = len(found_combinations)
    print(f"\nTotal number of combinations found: {count}")
    
    # Final answer in the specified format.
    print(f"<<<{count}>>>")

if __name__ == "__main__":
    find_fibonacci_sum_combinations()