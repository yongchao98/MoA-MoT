def solve_work_equation():
    """
    This function prints the derived formula for the work done by the current source.
    """
    # Define symbolic variables for clarity in the final expression
    mu = "μ"
    mu_0 = "μ₀"
    N = "N"
    w = "w"
    x1 = "x₁"
    x2 = "x₂"
    g = "g"
    I1 = "I₁"
    I2 = "I₂"

    # The final derived equation for the work done (W) by the current source
    # is formatted to be clear and readable.
    
    print("The work 'W' done by the current source for each cycle is given by the equation:")
    
    # The equation is constructed using an f-string for clear formatting.
    # The structure matches option D from the problem description.
    # Each part of the equation is explicitly shown.
    numerator_part1 = f"({mu} - {mu_0})"
    numerator_part2 = f"{N}²"
    numerator_part3 = f"{w}"
    numerator_part4 = f"({x2} - {x1})"
    numerator_part5 = f"({I2}² - {I1}²)"
    denominator = f"2{g}"

    print(f"W = - ( {numerator_part1} / {denominator} ) * {numerator_part2} * {numerator_part3} * {numerator_part4} * {numerator_part5}")

# Execute the function to print the result
solve_work_equation()