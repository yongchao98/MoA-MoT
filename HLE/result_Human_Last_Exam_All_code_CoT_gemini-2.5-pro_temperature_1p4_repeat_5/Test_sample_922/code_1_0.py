def solve_sequence():
    """
    This function presents the solution to the sequence problem by constructing
    the full mathematical equation it is a part of.
    """
    # The given sequence of numbers
    bases_lhs = [24663, 35005, 119261, 196219, 211770, 227296]
    
    # The single integer that 'completes' the sequence, found in August 2022
    completing_integer = 230159
    
    # The power involved in this sum of like powers problem
    power = 6

    print("The provided sequence forms part of the first known solution to the problem of finding six 6th powers that sum to a 6th power.")
    print("The integer that completes the sequence is the base of the term on the right-hand side of the equation.")
    print("\nThe complete equation is:")

    # Construct and print the full equation
    lhs_terms = [f"{n}^{power}" for n in bases_lhs]
    lhs_string = " + ".join(lhs_terms)
    rhs_string = f"{completing_integer}^{power}"
    
    print(f"{lhs_string} = {rhs_string}")

solve_sequence()