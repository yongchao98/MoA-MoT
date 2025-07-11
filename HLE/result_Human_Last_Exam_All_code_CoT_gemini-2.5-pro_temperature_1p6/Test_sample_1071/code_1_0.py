def find_bounds_of_t():
    """
    This function explains the derivation to find the lower and upper bounds of t.
    """
    
    print("Step 1: Define the problem in terms of intervals.")
    print("Let A be the interval [-1, t]. The variables a0, a1, a2, a3 are in A.")
    print("Let S be the set of all possible sums a_i + a_j, where a_i, a_j are in A.")
    print("S is the interval [-2, 2*t]. For A to be valid, t >= -1.")
    print("The requirement is: for any s0 in S, there must exist s1 in S such that s0 * s1 = 1.")
    print("-" * 50)

    print("Step 2: Handle the case where the sum is zero.")
    print("If s0 = 0, the equation 0 * s1 = 1 cannot be satisfied.")
    print("Therefore, 0 must not be in the interval S = [-2, 2*t].")
    print("This forces 2*t < 0, which means t < 0.")
    print("Combining with t >= -1, the valid range for t is [-1, 0).")
    print("-" * 50)

    print("Step 3: Analyze the condition for t < 0.")
    print("The condition means the set {1/s | s in S} must be a subset of S.")
    print("For t in [-1, 0), S = [-2, 2t].")
    print("The set of reciprocals, 1/S, is [1/(2t), 1/(-2)] = [1/(2t), -1/2].")
    print("We need [1/(2t), -1/2] to be a subset of [-2, 2t].")
    print("-" * 50)

    print("Step 4: Set up and solve the inequalities.")
    print("The subset condition gives two inequalities:")
    ineq1_str = "1/(2t) >= -2"
    ineq2_str = "-1/2 <= 2t"
    
    print(f"1) Lower bounds: {ineq1_str}")
    print("   Since t is negative, multiplying by 2t reverses the inequality sign:")
    print("   1 <= -2 * (2t)  =>  1 <= -4t  =>  t <= -1/4")
    
    print(f"\n2) Upper bounds: {ineq2_str}")
    print("   Dividing by 2 gives: t >= -1/4")
    print("-" * 50)

    print("Step 5: Conclusion.")
    print("The two conditions t <= -1/4 and t >= -1/4 can only be met if t = -1/4.")
    
    lower_bound = -0.25
    upper_bound = -0.25
    
    print("\nThe problem is solved at the boundary conditions, where the inequalities become equalities.")
    print("Final equations for the bounds, with t = -1/4:")
    
    # "output each number in the final equation"
    # Showing how t = -0.25 satisfies the boundary conditions as equalities.
    t_val = -0.25
    eq1_val1 = 2 * t_val
    eq1_val2 = -0.5
    print(f"For the upper bound: 2*t = {eq1_val1}, which must equal -1/2 = {eq1_val2}.")
    
    eq2_val1 = 1 / (2 * t_val)
    eq2_val2 = -2
    print(f"For the lower bound: 1/(2t) = {eq2_val1}, which must equal -2 = {eq2_val2}.")

    print("\nThe only possible value for t is -0.25.")
    print("Therefore, the lower and upper bounds are identical.")
    print(f"\n{lower_bound} {upper_bound}")

find_bounds_of_t()