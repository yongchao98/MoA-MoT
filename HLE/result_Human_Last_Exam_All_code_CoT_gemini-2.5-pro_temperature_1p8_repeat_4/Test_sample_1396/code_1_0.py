def solve_cake_cutting_bound():
    """
    Calculates the upper bound on the number of cuts for a connected
    epsilon-envy-free cake allocation for n agents using Stromquist's formula.
    """
    # The number of agents in the problem.
    n = 4

    # The formula for the upper bound (O) is n * (n - 1) / 2.
    # We will calculate this step by step.
    
    term_in_parentheses = n - 1
    numerator = n * term_in_parentheses
    denominator = 2
    
    # Perform the final calculation.
    upper_bound = numerator / denominator

    # Print the explanation and the step-by-step calculation.
    print(f"The upper bound O for a connected epsilon-envy-free allocation for n={n} agents is calculated using the formula: n * (n - 1) / 2")
    print("The calculation steps are:")
    print(f"O = ({n} * ({n} - 1)) / {denominator}")
    print(f"O = ({n} * {term_in_parentheses}) / {denominator}")
    print(f"O = {numerator} / {denominator}")
    print(f"O = {int(upper_bound)}")

solve_cake_cutting_bound()