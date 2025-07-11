def solve():
    """
    Calculates the sum of squares of coefficients for the given polynomial expansion.
    The problem is solved by setting up a recurrence relation which is then solved
    using two integer sequences x_n and z_n.
    The final answer is given by the formula 2 * x_19 + 5 * z_19.
    """
    
    n = 19
    
    # Initialize the sequences x_n and z_n based on their definitions.
    # x_n corresponds to lambda_2^n + lambda_3^n
    # z_n corresponds to (lambda_2^n - lambda_3^n) / sqrt(7)
    # where lambda_{2,3} = 3 +- sqrt(7)
    
    x_prev, x_curr = 2, 6  # Corresponds to x_0, x_1
    z_prev, z_curr = 0, 2  # Corresponds to z_0, z_1
    
    # Handle the base cases n=0 and n=1
    if n == 0:
        x_n = x_prev
        z_n = z_prev
    elif n == 1:
        x_n = x_curr
        z_n = z_curr
    else:
        # Loop to compute x_19 and z_19
        for _ in range(2, n + 1):
            x_next = 6 * x_curr - 2 * x_prev
            z_next = 6 * z_curr - 2 * z_prev
            
            x_prev, x_curr = x_curr, x_next
            z_prev, z_curr = z_curr, z_next
            
        x_n = x_curr
        z_n = z_curr

    # The final result is derived from the matrix calculation
    final_result = 2 * x_n + 5 * z_n
    
    print("The final result is calculated using the formula: 2 * x_19 + 5 * z_19.")
    print(f"The value for x_19 is: {x_n}")
    print(f"The value for z_19 is: {z_n}")
    print(f"The final equation is: 2 * {x_n} + 5 * {z_n} = {final_result}")

solve()