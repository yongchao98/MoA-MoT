def display_d_threshold_formula():
    """
    This function displays the formula for the d-threshold for Hamiltonicity.
    The formula is derived from known results in random graph theory.
    """

    # The formula for the d-threshold p is p = 1/2 - d/n.
    # We define the numbers and variables to construct and print the final equation as requested.
    
    number_one = 1
    number_two = 2
    variable_d = "d"
    variable_n = "n"

    print("The d-threshold for Hamiltonicity is given by the equation:")
    print(f"p = {number_one}/{number_two} - {variable_d}/{variable_n}")

if __name__ == '__main__':
    display_d_threshold_formula()