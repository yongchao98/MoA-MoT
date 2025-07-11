def solve():
    """
    Calculates the sum of the squares of the coefficients of the given polynomial expansion.
    The problem boils down to solving the recurrence relation x_n = 6*x_{n-1} - 2*x_{n-2}
    with initial conditions x_0 = 4 and x_1 = 22. We need to find x_19.
    """
    
    # We need to compute up to the 20th term (index 19).
    n = 20
    x = [0] * n
    
    # Initial conditions
    x[0] = 4
    x[1] = 22
    
    # Apply the recurrence relation
    for i in range(2, n):
        x[i] = 6 * x[i-1] - 2 * x[i-2]
        
    # The problem asks to output each number in the final equation.
    # The final equation is x[19] = 6 * x[18] - 2 * x[17].
    # The result is x[19].
    
    final_val = x[n-1]
    prev_val = x[n-2]
    prev_prev_val = x[n-3]
    
    print(f"The value is found by solving the recurrence x_n = 6*x_{n-1} - 2*x_{n-2}.")
    print(f"The initial values are x_0 = {x[0]} and x_1 = {x[1]}.")
    print(f"The final calculation is:")
    print(f"x_19 = 6 * x_18 - 2 * x_17")
    print(f"x_19 = 6 * {prev_val} - 2 * {prev_prev_val}")
    print(f"x_19 = {final_val}")

solve()