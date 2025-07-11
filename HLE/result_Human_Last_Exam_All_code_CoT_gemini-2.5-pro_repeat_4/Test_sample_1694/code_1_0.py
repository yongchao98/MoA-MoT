import cmath

def solve_and_sum_first_coordinate():
    """
    Solves the two independent matrix equations to find the sum of the
    first coordinate of the solutions.
    """

    # --- First Equation ---
    # The first equation is:
    # (diag(5, -5)) * X1^2 + X1^2 * (diag(6, 6)) = diag(-53/12, 0)
    # Let Y1 = X1^2. The equation for the (1,1) element y1_11 is:
    # 5 * y1_11 + y1_11 * 6 = -53/12
    # 11 * y1_11 = -53/12
    
    y1_11 = -53 / (12 * 11)
    
    # The (1,1) element of X1, let's call it x1_11, satisfies x1_11^2 = y1_11.
    # The solutions for x1_11 are the two square roots of y1_11.
    sol1_a = cmath.sqrt(y1_11)
    sol1_b = -cmath.sqrt(y1_11)
    
    # The sum of the first coordinates for the first equation's solutions
    sum1 = sol1_a + sol1_b
    
    # --- Second Equation ---
    # The second equation is:
    # (diag(4, -5)) * X2^2 + X2^2 * (diag(6, 6)) = diag(-3/11, 0)
    # Let Y2 = X2^2. The equation for the (1,1) element y2_11 is:
    # 4 * y2_11 + y2_11 * 6 = -3/11
    # 10 * y2_11 = -3/11
    
    y2_11 = -3 / (11 * 10)
    
    # The (1,1) element of X2, let's call it x2_11, satisfies x2_11^2 = y2_11.
    # The solutions for x2_11 are the two square roots of y2_11.
    sol2_a = cmath.sqrt(y2_11)
    sol2_b = -cmath.sqrt(y2_11)
    
    # The sum of the first coordinates for the second equation's solutions
    sum2 = sol2_a + sol2_b
    
    # --- Total Sum ---
    # The total sum is the sum of the sums from both equations.
    total_sum = sum1 + sum2
    
    # The sums will be 0j (complex zero). We take the real part for cleaner printing.
    sum1_real = sum1.real
    sum2_real = sum2.real
    total_sum_real = total_sum.real
    
    # Print the final equation as requested.
    # We show the sum for each equation and the final total sum.
    print(f"The sum of the first coordinates for the first equation is: {sum1_real}")
    print(f"The sum of the first coordinates for the second equation is: {sum2_real}")
    print("The final summation is:")
    print(f"{int(sum1_real)} + {int(sum2_real)} = {int(total_sum_real)}")

# Execute the function
solve_and_sum_first_coordinate()