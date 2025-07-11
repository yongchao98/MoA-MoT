def solve():
    """
    Finds the set of all fixed points and cycle values for the given number transformation process.
    """

    def transform(n):
        """
        Calculates B - A + 1 for a given number n.
        A is the smallest number formed from n's digits.
        B is the largest number formed from n's digits.
        """
        s = str(n)
        
        # Sort digits to find A and B
        sorted_digits = sorted(list(s))
        a_str = "".join(sorted_digits)
        b_str = "".join(sorted_digits[::-1])
        
        A = int(a_str)
        B = int(b_str)
        
        return B - A + 1

    # After one step, any 3-digit number transforms into one of these 10 values.
    # (d_max - d_min) can be 0 through 9.
    # N_new = 99 * (d_max - d_min) + 1
    seeds = {99 * i + 1 for i in range(10)}

    cycle_and_fixed_point_values = set()
    
    # We only need to trace the paths from these 10 seed values.
    for seed in seeds:
        path = []
        current_num = seed
        
        # Follow the path until we find a number we've already seen in this path.
        while current_num not in path:
            path.append(current_num)
            current_num = transform(current_num)
        
        # A cycle is found. Find where it starts and add all its members to our set.
        cycle_start_index = path.index(current_num)
        cycle = path[cycle_start_index:]
        cycle_and_fixed_point_values.update(cycle)
        
    # Sort the final set of values in ascending order for the output.
    sorted_values = sorted(list(cycle_and_fixed_point_values))
    
    # Format the output as a set string.
    result_str = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_str)

solve()