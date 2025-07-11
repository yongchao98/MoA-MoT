def solve_digit_process():
    """
    This script finds the total set of values that are either fixed points or part of a cycle
    when a specific iterative process is run on all positive three-digit numbers.

    The process is as follows:
    1. Take a number N.
    2. Form A, the smallest number from N's digits, and B, the largest.
    3. The new number is B - A + 1.
    4. Repeat with the new number.
    """

    def get_next_value(n):
        """
        Calculates the next value in the sequence for a given number n.
        For example, for 312, A=123, B=321, and the function returns B-A+1 = 199.
        """
        s_n = str(n)
        
        # Sort digits ascending to form number A
        a_digits = sorted(list(s_n))
        a = int("".join(a_digits))
        
        # Sort digits descending to form number B
        b_digits = a_digits[::-1]
        b = int("".join(b_digits))
        
        return b - a + 1

    # This set will store all unique numbers that are part of a cycle or are a fixed point.
    cycle_and_fixed_points = set()

    # Iterate through all positive three-digit numbers (100 to 999).
    for i in range(100, 1000):
        current_num = i
        path_list = []  # Tracks the sequence of numbers to identify the cycle later.
        path_set = set()  # For fast checking if a number has been seen in the current path.

        # Generate the sequence until a number repeats.
        while current_num not in path_set:
            path_list.append(current_num)
            path_set.add(current_num)
            current_num = get_next_value(current_num)
        
        # A cycle is detected. The cycle starts at the first occurrence of current_num.
        try:
            cycle_start_index = path_list.index(current_num)
            # The cycle is the sublist from the start index to the end.
            cycle = path_list[cycle_start_index:]
            # Add all numbers in this cycle to our main set.
            cycle_and_fixed_points.update(cycle)
        except ValueError:
            # This case is not expected due to the logic of the while loop.
            pass

    # Sort the final set of numbers in ascending order.
    sorted_values = sorted(list(cycle_and_fixed_points))

    # Format the output string as requested, e.g., {1, 2, 3}.
    result_str = "{" + ", ".join(map(str, sorted_values)) + "}"
    print(result_str)

solve_digit_process()