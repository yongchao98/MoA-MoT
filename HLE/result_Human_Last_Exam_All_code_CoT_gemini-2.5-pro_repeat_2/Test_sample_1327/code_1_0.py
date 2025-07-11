def find_next_number():
    """
    This function finds the next number in the sequence based on the identified pattern.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    
    # The pattern is x_n = 3 * x_{n-1} - d_n
    # The sequence of subtracted numbers (d_n) is 10, 18, 34.
    d_sequence = [10, 18, 34]
    
    # The differences in the d_sequence are 8 (18-10) and 16 (34-18).
    # This difference doubles each time.
    last_diff = d_sequence[2] - d_sequence[1]  # 16
    next_diff = last_diff * 2  # 32
    
    # The next number to subtract is the last one plus the new difference.
    next_d = d_sequence[-1] + next_diff  # 34 + 32 = 66
    
    # The last known number in the original sequence.
    last_term = sequence[-1] # 119
    
    # Calculate the next number in the sequence.
    next_term = 3 * last_term - next_d # 3 * 119 - 66
    
    print("The pattern is that the next number is 3 times the previous number, minus a value.")
    print("The sequence of subtracted values is 10, 18, 34, ...")
    print("The differences between these values are 8, 16, ... which are doubling.")
    print("The next difference is 16 * 2 = 32.")
    print("So the next value to subtract is 34 + 32 = 66.")
    print(f"The final equation is: {3} * {119} - {66} = {next_term}")

find_next_number()