def solve_sequence():
    """
    Solves the number sequence puzzle and explains the pattern.
    """
    x = [2, 11, 23, 51, 119]
    
    # The pattern identified is x_{n+1} = 3*x_n - k_n starting from n=2.
    # Let's find the sequence k_n.
    k2 = 3 * x[1] - x[2]  # k for x_3
    k3 = 3 * x[2] - x[3]  # k for x_4
    k4 = 3 * x[3] - x[4]  # k for x_5
    
    # The sequence of subtrahends is k = [10, 18, 34, ...]
    
    # The differences in k are:
    diff1 = k3 - k2  # 18 - 10 = 8
    diff2 = k4 - k3  # 34 - 18 = 16
    
    # The differences form a geometric progression with ratio 2.
    # The next difference is:
    next_diff = diff2 * 2  # 16 * 2 = 32
    
    # The next subtrahend is:
    next_k = k4 + next_diff # 34 + 32 = 66
    
    # Now, calculate the next term in the main sequence.
    next_x = 3 * x[4] - next_k
    
    print("The sequence is: 2, 11, 23, 51, 119, ...")
    print("The pattern is that each term (starting from the third) is 3 times the previous term minus a value.")
    print(f"Equation for the 3rd term: {x[2]} = 3 * {x[1]} - {k2}")
    print(f"Equation for the 4th term: {x[3]} = 3 * {x[2]} - {k3}")
    print(f"Equation for the 5th term: {x[4]} = 3 * {x[3]} - {k4}")
    print("\nThe sequence of subtrahends is 10, 18, 34, ...")
    print("The differences in this sequence are 8, 16, ... which double each time.")
    print(f"The next difference is {diff2} * 2 = {next_diff}.")
    print(f"The next subtrahend is {k4} + {next_diff} = {next_k}.")
    print("\nTherefore, the final equation to find the missing number is:")
    print(f"Missing Number = 3 * {x[4]} - {next_k} = {next_x}")

solve_sequence()
<<<A>>>