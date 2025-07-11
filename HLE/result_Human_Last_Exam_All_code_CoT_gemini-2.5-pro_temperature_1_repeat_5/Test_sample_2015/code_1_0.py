def solve_cycle_problem():
    """
    This function calculates the set of all fixed values and cycle members
    for the described process, starting with all three-digit numbers.
    """

    def calculate_next_value(n):
        """
        Calculates the B - A + 1 value for a given number n.
        A is the smallest number formed by n's digits.
        B is the largest number formed by n's digits.
        """
        s_n = str(n)
        
        # Form A, the smallest number
        # Note: int() correctly handles leading zeros, e.g., int("012") -> 12
        a_str = "".join(sorted(list(s_n)))
        a = int(a_str)

        # Form B, the largest number
        b_str = "".join(sorted(list(s_n), reverse=True))
        b = int(b_str)

        return b - a + 1

    # This set will store all numbers that are part of a cycle or are a fixed point.
    final_cycle_values = set()

    # Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        history = []
        current_num = start_num

        # Generate the sequence until a number is repeated.
        while current_num not in history:
            history.append(current_num)
            current_num = calculate_next_value(current_num)

        # A cycle is found. Identify the members of the cycle.
        cycle_start_index = history.index(current_num)
        cycle = history[cycle_start_index:]
        
        # Add the cycle members to our final set.
        final_cycle_values.update(cycle)

    # Sort the final values and print in the required format.
    sorted_values = sorted(list(final_cycle_values))
    
    # The prompt asks for outputting the numbers in the final result.
    # The following print statement formats the result as a set string.
    result_str = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_str)

# Run the solver function.
solve_cycle_problem()