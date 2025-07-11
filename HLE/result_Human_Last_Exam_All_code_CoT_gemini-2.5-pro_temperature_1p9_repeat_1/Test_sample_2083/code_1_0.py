def solve():
    """
    This function determines and prints the minimum hidden-layer width
    required to compute the squared norm of an N-dimensional vector.
    """

    # The dimension of the input vector is denoted by N.
    # We treat it as a symbolic character for the output formula.
    variable_name = 'N'

    # The analysis shows the minimum required width is 2 times N.
    # The number in this equation is 2.
    coefficient = 2
    
    # We print the final formula. The problem asks to output the number(s) in the equation.
    # In '2 * N', the number is 2.
    print(f"{coefficient} * {variable_name}")

solve()