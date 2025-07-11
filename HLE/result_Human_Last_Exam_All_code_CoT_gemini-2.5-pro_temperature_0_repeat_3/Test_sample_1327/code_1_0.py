def solve_sequence():
    """
    This function solves the given number sequence puzzle.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    
    # The pattern is a_n = 3 * a_{n-1} - b_n for n >= 3
    # Let's find the sequence b_n
    # b_1 corresponds to a_3, b_2 to a_4, etc.
    b = []
    for i in range(2, len(sequence)):
        b_i = 3 * sequence[i-1] - sequence[i]
        b.append(b_i)
        
    # b is [10, 18, 34]
    
    # Now we find the pattern in b. The difference between terms doubles.
    # diff1 = 18 - 10 = 8
    # diff2 = 34 - 18 = 16
    # The next difference will be 16 * 2 = 32
    last_diff = b[-1] - b[-2]
    next_diff = last_diff * 2
    
    # The next term in b is the last term of b + the next difference
    next_b = b[-1] + next_diff
    
    # Now we can calculate the next term in the main sequence
    last_term_in_sequence = sequence[-1]
    next_term = 3 * last_term_in_sequence - next_b
    
    print("The sequence is: 2, 11, 23, 51, 119, ...")
    print("The pattern is a_n = 3 * a_{n-1} - b_n, where b_n is a sequence starting from the 3rd term.")
    print(f"The sequence of subtracted numbers 'b' is: {b}")
    print("The pattern in 'b' is that the difference between consecutive terms doubles: 8, 16, ...")
    print(f"The next difference is 16 * 2 = {next_diff}")
    print(f"The next number to subtract is {b[-1]} + {next_diff} = {next_b}")
    print("\nTherefore, the final equation is:")
    print(f"3 * {last_term_in_sequence} - {next_b} = {next_term}")

solve_sequence()