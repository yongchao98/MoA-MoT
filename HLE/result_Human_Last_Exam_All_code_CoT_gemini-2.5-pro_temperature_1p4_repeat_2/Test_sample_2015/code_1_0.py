def find_three_digit_cycles():
    """
    This function calculates the set of all values that are part of a fixed point
    or a cycle for the process B-A+1 on all three-digit numbers.
    """
    N = 3  # The number of digits for the process
    all_cycle_members = set()

    # 1. Iterate through all positive three-digit numbers.
    for start_num in range(100, 1000):
        history = []
        current_num = start_num

        # 2. Simulate the process until a number is repeated.
        while current_num not in history:
            history.append(current_num)

            # Ensure the number is treated as having N digits, padding with zeros if needed.
            s_num = str(current_num).zfill(N)

            # Create A (smallest) and B (largest) from the digits.
            sorted_digits = sorted(list(s_num))
            
            a_str = "".join(sorted_digits)
            b_str = "".join(sorted_digits[::-1]) # Reverse for descending order

            a = int(a_str)
            b = int(b_str)

            # Calculate the next number in the sequence.
            current_num = b - a + 1
        
        # 3. A cycle is found. Identify its members.
        # The cycle starts at the first index of the repeated number.
        try:
            cycle_start_index = history.index(current_num)
            cycle = history[cycle_start_index:]
            
            # 4. Add the cycle members to the master set.
            all_cycle_members.update(cycle)
        except ValueError:
            # This case should not be reached due to the loop's logic.
            pass

    # 5. Sort the final collection of numbers.
    sorted_values = sorted(list(all_cycle_members))

    # Print the final result in the specified format.
    print(f"{{{', '.join(map(str, sorted_values))}}}")

if __name__ == '__main__':
    find_three_digit_cycles()