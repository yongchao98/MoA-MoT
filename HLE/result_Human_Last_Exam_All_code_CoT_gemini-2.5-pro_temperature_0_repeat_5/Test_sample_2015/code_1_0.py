def solve_digit_reordering_cycle():
    """
    This function simulates a digit reordering process for all three-digit numbers
    to find the total set of values that are part of terminal cycles or fixed points.

    The process for a number n is:
    1. Form A, the smallest number from n's digits, and B, the largest.
    2. The next number is B - A + 1.
    3. This is repeated until a value repeats.
    """
    # memo will store the cycle/fixed point set for any number already processed.
    memo = {}
    final_set = set()

    # 1. Iterate through all positive three-digit numbers.
    for i in range(100, 1000):
        # If we have already found the cycle for this number's path, we can skip it.
        if i in memo:
            continue

        path = []
        path_set = set()
        current_num = i

        # 2. Trace the path until we find a number we've seen before.
        while current_num not in memo and current_num not in path_set:
            path.append(current_num)
            path_set.add(current_num)

            s_num = str(current_num)
            
            # Form A (smallest) and B (largest) from the digits.
            digits = sorted(list(s_num))
            A = int("".join(digits))
            B = int("".join(reversed(digits)))
            
            # Calculate the next number in the sequence.
            current_num = B - A + 1

        # 3. Determine the cycle/fixed point.
        if current_num in memo:
            # The path leads to a known cycle.
            cycle = memo[current_num]
        else:
            # A new cycle was found within the current path.
            cycle_start_index = path.index(current_num)
            cycle = set(path[cycle_start_index:])

        # 4. Add the numbers from the found cycle to our final set.
        final_set.update(cycle)

        # Memoize the result for all numbers in the path we just traversed.
        for num in path:
            memo[num] = cycle

    # 5. Sort the final set of numbers in ascending order for the output.
    sorted_results = sorted(list(final_set))
    
    # Print the result in the required format: {val1, val2, ...}
    print(f"{{{', '.join(map(str, sorted_results))}}}")

solve_digit_reordering_cycle()