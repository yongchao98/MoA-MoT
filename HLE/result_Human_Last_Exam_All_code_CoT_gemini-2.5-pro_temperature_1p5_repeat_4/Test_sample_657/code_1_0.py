def solve_controller():
    """
    This function prints the set of all proper stabilizing controllers H_2(s)
    for the given plant H_1(s).
    The result is presented as a transfer function parametrized by K(s),
    a free parameter which can be any stable and proper rational function.
    """

    # The numerator is: 4*s^2 + 8*s + 4 + (s^2 - 1)*K(s)
    # The denominator is: s^2 - 1 - s*K(s)

    num_part1 = "(4)*s^2 + (8)*s + (4)"
    num_part2 = "((1)*s^2 + (-1))*K(s)"
    
    den_part1 = "(1)*s^2 + (-1)"
    den_part2 = "(1)*s*K(s)"
    
    # Constructing the string representation for clarity
    line_length = len(f"H_2(s) =  {num_part1} + {num_part2} ")
    separator = "-" * line_length
    
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("")
    print(f"         {num_part1} + {num_part2}")
    print(f"H_2(s) = {separator}")
    print(f"         {den_part1} - {den_part2}")
    print("")
    print("where K(s) is any stable and proper rational function.")

if __name__ == "__main__":
    solve_controller()