def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a multi-layered pattern.
    1. It establishes a primary recurrence relation: a_{n+1} = 2*a_n + c_n.
    2. It calculates the sequence of constants, c_n.
    3. It finds a pattern within the c_n sequence by examining its differences, d_n.
    4. It extrapolates the next terms for d_n, c_n, and finally the original sequence a_n.
    """
    # The given sequence
    a = [2, 11, 23, 51, 119]
    
    # Let's analyze the pattern a_{n+1} = 2*a_n + c_n, focusing on n >= 2
    c2 = a[2] - 2 * a[1]  # 23 - 2*11 = 1
    c3 = a[3] - 2 * a[2]  # 51 - 2*23 = 5
    c4 = a[4] - 2 * a[3]  # 119 - 2*51 = 17
    
    # The sequence of constants (starting from the second term) is 1, 5, 17
    # Let's find the pattern in this new sequence by looking at the differences
    d2 = c3 - c2  # 5 - 1 = 4
    d3 = c4 - c3  # 17 - 5 = 12
    
    # The differences {4, 12} form a geometric progression with ratio 3
    ratio = d3 // d2
    
    # The next difference d4 follows this pattern
    d4 = d3 * ratio  # 12 * 3 = 36
    
    # Now, we find the next constant c5
    c5 = c4 + d4  # 17 + 36 = 53
    
    # Finally, we find the next term in the original sequence, a6
    a5 = a[4]
    a6 = 2 * a5 + c5  # 2 * 119 + 53 = 291
    
    # The final equation and its result
    # We output each number as requested
    num1 = 2
    num2 = a5
    num3 = c5
    result = a6
    
    print(f"{num1} * {num2} + {num3} = {result}")

solve_sequence()