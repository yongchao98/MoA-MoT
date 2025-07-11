import sympy

def solve_bounds():
    """
    This function explains the logical steps to find the lower and upper bounds for t.
    """
    print("Step 1: Define the problem and variables.")
    # The variables a0, a1, a2, a3 are in the interval [-1, t].
    # The requirement is: (a0 + a2) * (a1 + a3) = 1.
    # The number in this equation is 1.
    print("The governing equation is (a_0 + a_2)(a_1 + a_3) = 1.")
    
    print("\nStep 2: Determine the range of the sums.")
    # Let S = a_i + a_j. The interval for a_i is [-1, t].
    # For this interval to be non-empty, we must have t >= -1.
    # The minimum sum is -1 + (-1) = -2. The maximum sum is t + t = 2*t.
    # So the range of possible sums is R = [-2, 2*t].
    # The numbers defining the basic sum range are -2 and 2.
    print("The sum of any two numbers from [-1, t] lies in the interval [-2, 2*t].")

    print("\nStep 3: Analyze the condition S=0.")
    # Let S1 = a_0 + a_2. The equation requires finding a_1, a_3 such that a_1 + a_3 = 1/S1.
    # If S1 can be 0, 1/S1 is undefined and no solution exists.
    # S1 can be 0 if the range [-2, 2*t] contains 0. This happens if 2*t >= 0, so t >= 0.
    # If t >= 0, we can always choose a_0 and a_2 to make their sum 0 (e.g., a_0=x, a_2=-x).
    # Therefore, t must be strictly less than 0 to avoid S1=0.
    print("To avoid division by zero, the sum a_0 + a_2 cannot be 0.")
    print("This implies t must be less than 0.")

    print("\nStep 4: Solve for t when t < 0.")
    # We now have -1 <= t < 0.
    # The range of sums R = [-2, 2*t] contains only negative numbers.
    # The condition is: for any s in R, 1/s must also be in R.
    # Let's find the image of R under the function f(s) = 1/s.
    # f([-2, 2*t]) = [1/(2*t), 1/(-2)] = [1/(2*t), -0.5].
    # This image must be a subset of R: [1/(2*t), -0.5] subset of [-2, 2*t].
    print("For t < 0, we require that for any s in [-2, 2t], 1/s is also in [-2, 2t].")
    print("This leads to two inequalities based on the interval boundaries.")
    
    t = sympy.Symbol('t')
    
    # Inequality 1: Lower bound
    # 1/(2*t) >= -2
    # Since t < 0, 2t is negative, so we flip the inequality sign on multiplying.
    # 1 <= -4*t  =>  t <= -1/4
    ineq1_lhs = 1/(2*t)
    ineq1_rhs = -2
    print(f"Inequality 1: {ineq1_lhs} >= {ineq1_rhs}")
    # We solve for t, knowing t is negative.
    solution1 = sympy.solve(1 <= -4*t, t)
    print(f"Solving for t gives: {solution1}")

    # Inequality 2: Upper bound
    # -0.5 <= 2*t
    # t >= -0.5 / 2 => t >= -0.25
    ineq2_lhs = -0.5
    ineq2_rhs = 2*t
    print(f"Inequality 2: {ineq2_lhs} <= {ineq2_rhs}")
    solution2 = sympy.solve(ineq2_lhs <= ineq2_rhs, t)
    print(f"Solving for t gives: {solution2}")
    
    # The numbers in these final inequalities are -1 and 4.

    print("\nStep 5: Combine results and conclude.")
    # Combining t <= -1/4 and t >= -1/4 gives t = -1/4.
    lower_bound = -0.25
    upper_bound = -0.25
    print(f"The only value of t satisfying both inequalities is t = {lower_bound}.")
    print("Therefore, the lower and upper bounds for t are the same.")
    
    print("\nFinal Answer:")
    print(f"The lower and upper bounds of t are: {lower_bound} {upper_bound}")

solve_bounds()
<<<-0.25 -0.25>>>