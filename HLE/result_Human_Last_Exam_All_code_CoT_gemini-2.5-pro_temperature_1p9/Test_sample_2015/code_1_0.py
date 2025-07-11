import sys

def solve_and_print():
    """
    Solves the problem of finding all cycle and fixed-point values for a given process
    on three-digit numbers, and prints the result in the specified format.
    """

    def calculate_next_value(n):
        """
        Calculates the next value in the sequence for a given number n.
        It computes B (largest number from digits of n) and A (smallest number),
        and returns B - A + 1.
        """
        s = str(n)
        s_sorted_list = sorted(s)
        
        # Form the smallest number, A
        a_str = "".join(s_sorted_list)
        a = int(a_str)
        
        # Form the largest number, B
        # The list is already sorted, reversing it gives the largest
        b_str = "".join(reversed(s_sorted_list))
        b = int(b_str)
        
        return b - a + 1

    # This set will store all numbers that are part of a cycle or are fixed points.
    global_cycle_members = set()

    # Iterate through all positive three-digit numbers.
    for start_n in range(100, 1000):
        
        path = []
        seen_in_path = set()
        current_n = start_n
        
        # Trace the sequence starting from start_n until a value repeats.
        while current_n not in seen_in_path:
            seen_in_path.add(current_n)
            path.append(current_n)
            current_n = calculate_next_value(current_n)
            
        # A cycle is found. The repeated number is `current_n`.
        # The cycle begins at the first occurrence of `current_n` in the path.
        try:
            cycle_start_index = path.index(current_n)
            cycle = path[cycle_start_index:]
            # Add all members of this cycle to our global set.
            global_cycle_members.update(cycle)
        except ValueError:
            # This should not be reached given the logic of the while loop.
            pass

    # Sort the final set of numbers in ascending order.
    sorted_values = sorted(list(global_cycle_members))

    # Format the result as a string in set notation.
    result_string = "{" + ", ".join(map(str, sorted_values)) + "}"
    
    print(result_string)

solve_and_print()