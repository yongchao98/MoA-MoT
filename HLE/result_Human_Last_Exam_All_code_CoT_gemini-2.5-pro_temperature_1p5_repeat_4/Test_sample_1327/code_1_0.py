def solve_sequence():
    """
    This function solves the number sequence puzzle.
    The sequence is 2, 11, 23, 51, 119, ...
    """
    a = [2, 11, 23, 51, 119]
    
    # The pattern is a_n = 3*a_{n-1} - k_{n-1} for n>2
    # The sequence k follows its own pattern: k_n = 2*k_{n-1} - 2
    # Let's define the initial k values derived from the main sequence
    # k_3 = (3*a_2) - a_3 = 3*11 - 23 = 10
    # k_4 = (3*a_3) - a_4 = 3*23 - 51 = 18
    # k_5 = (3*a_4) - a_5 = 3*51 - 119 = 34
    
    k = [10, 18, 34] # k values corresponding to a_3, a_4, a_5
    
    # Calculate the next k value
    next_k = 2 * k[-1] - 2
    
    # Calculate the next a value
    next_a = 3 * a[-1] - next_k
    
    # We still need to output the original numbers in the final equation.
    # The final equation is a_6 = 3 * a_5 - k_6
    print(f"The pattern is a_n = 3 * a_{{n-1}} - k_n, where k_n = 2 * k_{{n-1}} - 2.")
    print(f"The last number in the sequence is a_5 = {a[-1]}.")
    print(f"The last subtraction amount is k_5 = {k[-1]}.")
    print(f"The next subtraction amount is k_6 = 2 * {k[-1]} - 2 = {next_k}.")
    print(f"The next number in the sequence is a_6 = 3 * {a[-1]} - {next_k} = {next_a}.")
    
solve_sequence()
<<<A>>>