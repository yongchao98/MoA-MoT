import itertools

def solve_fibonacci_combinations():
    """
    Finds and prints all combinations of 3 to 7 Fibonacci numbers whose sum is also a Fibonacci number.
    The combinations are drawn from the first 40 Fibonacci numbers.
    """
    # Step 1: Generate Fibonacci numbers for combinations (first 40)
    fibs_for_combinations = [1, 1]
    for _ in range(38):
        fibs_for_combinations.append(fibs_for_combinations[-1] + fibs_for_combinations[-2])

    # Step 2: Generate a larger set of Fibonacci numbers for checking the sum
    # The max sum will be < F_45, so we generate up to F_45 for the lookup set.
    fibs_for_lookup = [1, 1]
    for _ in range(43):
        fibs_for_lookup.append(fibs_for_lookup[-1] + fibs_for_lookup[-2])
    fib_lookup_set = set(fibs_for_lookup)

    # Step 3: Use a set to store unique equations found
    found_equations = set()

    # Step 4: Iterate through combination lengths (3 to 7)
    for k in range(3, 8):
        # Generate combinations of size k
        for combo in itertools.combinations(fibs_for_combinations, k):
            current_sum = sum(combo)

            # Check if the sum is a Fibonacci number
            if current_sum in fib_lookup_set:
                # Create a canonical representation of the equation and add to the set
                sorted_combo = tuple(sorted(combo))
                equation = " + ".join(map(str, sorted_combo)) + f" = {current_sum}"
                found_equations.add(equation)

    # Step 5: Print the results in a sorted order for readability
    # Sort by the sum, then by the equation string itself
    sorted_results = sorted(
        list(found_equations),
        key=lambda x: (int(x.split('=')[1].strip()), x)
    )

    for eq in sorted_results:
        print(eq)

    count = len(found_equations)
    print(f"\nTotal number of combinations found: {count}")
    
    # Final answer in the required format
    print(f"<<<{count}>>>")

solve_fibonacci_combinations()