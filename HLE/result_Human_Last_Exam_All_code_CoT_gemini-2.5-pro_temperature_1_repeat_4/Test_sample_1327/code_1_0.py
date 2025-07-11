def solve_sequence():
    """
    This function solves the number sequence problem by identifying the pattern
    and calculating the next term.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]
    
    # The pattern is a(n+1) = 3 * a(n) - C, starting from the second term.
    # Let's find the sequence of C.
    c1 = 3 * sequence[1] - sequence[2]  # 3*11 - 23 = 10
    c2 = 3 * sequence[2] - sequence[3]  # 3*23 - 51 = 18
    c3 = 3 * sequence[3] - sequence[4]  # 3*51 - 119 = 34
    
    # The sequence of subtractions is [10, 18, 34].
    # Let's find the pattern in this new sequence.
    diff1 = c2 - c1  # 18 - 10 = 8
    diff2 = c3 - c2  # 34 - 18 = 16
    
    # The differences are doubling (8, 16, ...). The next difference is 32.
    next_diff = diff2 * 2
    
    # The next number in the subtraction sequence.
    next_c = c3 + next_diff
    
    # The last known term in the original sequence.
    last_term = sequence[-1]
    
    # Calculate the next term in the original sequence.
    next_term = 3 * last_term - next_c
    
    print("The pattern is that the next number is 3 times the current number minus a value.")
    print("This value is part of a sequence (10, 18, 34, ...) where the difference between terms doubles each time (8, 16, 32, ...).")
    print(f"The next value to subtract is {c3} + {next_diff} = {next_c}.")
    print("The final equation is:")
    print(f"3 * {last_term} - {next_c} = {next_term}")

solve_sequence()