import sympy

def solve_pde_lower_bound():
    """
    This function determines the lower bound 'a' for the given mathematical expression
    by analyzing a special case (a sharp interface solution) and prints the derivation.
    """

    # Define symbols for clarity in the explanation
    u = sympy.Symbol('u')
    u_bar = sympy.Symbol('ubar')
    u_bar_x = sympy.Symbol('ubar_x')
    x = sympy.Symbol('x')
    t = sympy.Symbol('t')
    a = sympy.Symbol('a')

    # Step 1: State the expression to be bounded
    print("Step 1: The expression to find a lower bound for is E.")
    print("E = (d/dt + (1-2u) * ubar * d/dx) * ubar\n")

    # Step 2: Consider a stationary sharp interface solution (step function)
    print("Step 2: We analyze a stationary sharp interface solution (a step function).")
    print("Let u(t,x) be a stationary step function: u(x) = 1 for x < 0, and u(x) = 0 for x > 0.")
    print("For a stationary solution, d/dt(u) = 0, which implies d/dt(ubar) = 0.\n")

    # Step 3: Simplify the expression for E
    print("Step 3: The expression for E simplifies to:")
    print("E = (1 - 2u) * ubar * ubar_x\n")

    # Step 4: Calculate the terms for this solution
    print("Step 4: Calculate the terms for this solution.")
    print("The derivative of the step function is the Dirac delta: u_x = -delta(x).")
    print("ubar is the convolution of the kernel (1/2)exp(-|x|) with u_x.")
    print("ubar(x) = - (1/2) * exp(-|x|)")
    print("The derivative of ubar is: ubar_x(x) = (1/2) * sgn(x) * exp(-|x|)\n")

    # Step 5: Evaluate E for x > 0 and x < 0
    print("Step 5: Evaluate E for the two regions x > 0 and x < 0.")
    print("Case 1: x > 0")
    print("u = 0, sgn(x) = 1")
    print("E(x) = (1 - 2*0) * (-1/2 * exp(-x)) * (1/2 * exp(-x))")
    print("E(x) = -1/4 * exp(-2x)\n")

    print("Case 2: x < 0")
    print("u = 1, sgn(x) = -1")
    print("E(x) = (1 - 2*1) * (-1/2 * exp(x)) * (-1/2 * exp(x))")
    print("E(x) = -1 * (1/4 * exp(2x))")
    print("E(x) = -1/4 * exp(2x)\n")

    # Step 6: Combine cases and find the infimum
    print("Step 6: Combine the results.")
    print("For any x != 0, the expression is E(x) = -1/4 * exp(-2|x|).")
    print("The infimum (greatest lower bound) of E(x) is the limit as x -> 0.\n")

    # Step 7: State the final answer
    lower_bound = -1/4
    print("Step 7: The lower bound 'a' is determined.")
    final_equation = sympy.Eq(a, lower_bound)
    
    # Print the final equation with each number
    a_val, bound_val = final_equation.args
    numerator, denominator = bound_val.as_numer_denom()
    print(f"The final equation for the lower bound is:")
    print(f"{a_val} = {numerator}/{denominator}")


if __name__ == "__main__":
    solve_pde_lower_bound()
    # The final answer in the required format
    # The value is -1/4 = -0.25
    print("\n<<<" + str(-1/4) + ">>>")