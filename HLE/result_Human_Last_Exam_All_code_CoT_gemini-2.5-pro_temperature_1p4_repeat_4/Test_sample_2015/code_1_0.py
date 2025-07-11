def find_cycle_values():
    """
    This script finds the total set of values included in fixed points or cycles
    for a specific iterative process applied to all three-digit numbers.

    The process is as follows:
    1. Start with a number N.
    2. Form A, the smallest number from N's digits, and B, the largest.
    3. The next number is B - A + 1.
    4. Repeat with the new number.

    The script iterates through all starting numbers from 100 to 999,
    detects the cycles or fixed points, and collects all unique values
    found within them.
    """
    
    # This set will store all numbers that are part of a fixed point or a cycle.
    all_cycle_members = set()

    # We iterate through all possible starting three-digit numbers.
    for start_n in range(100, 1000):
        
        n = start_n
        
        # history_list tracks the sequence of numbers to identify the cycle members.
        history_list = []
        # history_set provides a fast way to check if a number has been seen before.
        history_set = set()

        # The loop continues until we encounter a number that has been seen before.
        while n not in history_set:
            history_list.append(n)
            history_set.add(n)
            
            # Convert the number to a string to access its digits.
            s_n = str(n)
            
            # To get A (the smallest number), we sort the digits asc.
            s_a = "".join(sorted(s_n))
            a = int(s_a)
            
            # To get B (the largest number), we sort the digits desc.
            s_b = "".join(sorted(s_n, reverse=True))
            b = int(s_b)
            
            # Calculate the next number in the sequence.
            n = b - a + 1
        
        # 'n' is the first repeated number. Find where the cycle starts.
        try:
            cycle_start_index = history_list.index(n)
            # The cycle consists of all numbers from that index onwards.
            cycle = history_list[cycle_start_index:]
            # Add all members of this newly found cycle to our master set.
            all_cycle_members.update(cycle)
        except ValueError:
            # This case should not be reached with the current logic.
            # It's here for robustness.
            pass

    # Sort the final set of numbers in ascending order for the output.
    sorted_members = sorted(list(all_cycle_members))
    
    # Format the output as a set string, e.g., {1, 2, 3}.
    result_str = f"{{{', '.join(map(str, sorted_members))}}}"
    
    print(result_str)

find_cycle_values()
<<< {1, 100, 397, 496, 595} >>>