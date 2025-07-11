def find_cycle_values():
    """
    Finds the total set of values included in fixed points or cycles for the given process
    when run on all positive three-digit numbers.
    """

    # Memoization cache to store the determined cycle for any given number.
    # This avoids re-computing paths we have already explored.
    memo = {}

    # The final set of all numbers that are part of a cycle or fixed point.
    all_cycle_members = set()

    def transform(n):
        """
        Takes an integer n, treats it as a 3-digit number, and returns B - A + 1.
        """
        # Pad with leading zeros to ensure 3 digits. The result of the operation
        # is always <= 892, so it will not exceed 3 digits after the first step.
        s = str(n).zfill(3)

        # Sort the digits to get the smallest number A and largest number B.
        sorted_digits = sorted(s)
        a = int("".join(sorted_digits))
        b = int("".join(sorted_digits[::-1]))

        return b - a + 1

    # Iterate through every positive three-digit number.
    for start_num in range(100, 1000):
        # If we have already found the cycle for this number, skip it.
        if start_num in memo:
            continue

        # Track the sequence of numbers generated from the start_num.
        trajectory = []
        # A set for fast lookups to detect when a cycle begins.
        seen_in_trajectory = set()
        
        current_num = start_num

        # Loop until we find a cycle or hit a number we've already processed.
        while True:
            # If we hit a number for which the cycle is already known,
            # we can stop and use the cached result.
            if current_num in memo:
                cycle = memo[current_num]
                break

            # If we encounter a number that's already in the current path,
            # we have found a cycle.
            if current_num in seen_in_trajectory:
                cycle_start_index = trajectory.index(current_num)
                # The cycle is the part of the trajectory from the first encounter.
                cycle = set(trajectory[cycle_start_index:])
                break

            # If no cycle is found yet, continue the process.
            seen_in_trajectory.add(current_num)
            trajectory.append(current_num)
            current_num = transform(current_num)

        # Add the members of the found cycle to our master set.
        all_cycle_members.update(cycle)

        # Cache the result for every number in the path we just took.
        # This is the memoization step that speeds up the overall process.
        for n in trajectory:
            memo[n] = cycle

    # Sort the final set of numbers in ascending order for the output.
    sorted_values = sorted(list(all_cycle_members))
    
    # Print the result in the specified format.
    result_str = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_str)

find_cycle_values()
<<<
{100, 397, 496, 595}
>>>