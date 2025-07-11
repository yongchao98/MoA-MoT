def find_bounds_of_t():
    """
    This function determines the lower and upper bounds for the real number t
    based on the given mathematical requirement. It prints the derivation
    steps and the final answer.
    """
    print("Step 1: Define the set of sums.")
    print("Let a_i be in the interval [-1, t].")
    print("The sum S = a_i + a_j must lie in the interval [S_min, S_max].")
    s_min = -1 + -1
    print(f"S_min = -1 + (-1) = {s_min}")
    print("S_max = t + t = 2*t")
    print("The set of possible sums is I = [-2, 2*t].\n")

    print("Step 2: Analyze the condition.")
    print("The condition (a_0 + a_2)(a_1 + a_3) = 1 implies that for any sum s in I,")
    print("its reciprocal 1/s must also be in I.\n")

    print("Step 3: Exclude the possibility of division by zero.")
    print("If 0 is in I, we could have a sum of 0, for which 1/0 is undefined.")
    print("0 is in [-2, 2*t] if 2*t >= 0, which means t >= 0.")
    print("To prevent this, we must have t < 0.\n")

    print("Step 4: Set up inequalities for t < 0.")
    print("For t < 0, the set of reciprocals of I = [-2, 2*t] is J = [1/(2*t), -0.5].")
    print("The condition is that J must be a subset of I, i.e., [1/(2*t), -0.5] is a subset of [-2, 2*t].")
    print("This leads to two inequalities:")
    print("  1) Lower bound: -2 <= 1/(2*t)")
    print("  2) Upper bound: -0.5 <= 2*t\n")

    print("Step 5: Solve the inequalities.")
    # From -0.5 <= 2*t
    lower_bound = -0.5 / 2
    print(f"From -0.5 <= 2*t, we get t >= -0.5 / 2, so t >= {lower_bound}.")
    # From -2 <= 1/(2*t) => -4*t >= 1 (since t<0) => t <= -1/4
    upper_bound = -1.0 / 4
    print(f"From -2 <= 1/(2*t), we get -4*t >= 1, so t <= -1 / 4, so t <= {upper_bound}.")
    print("\nFor both inequalities to hold, t must be exactly -0.25.\n")
    
    print("Step 6: Conclusion and Example.")
    print(f"The only value for t is {lower_bound}.")
    print("Therefore, the lower and upper bounds are the same.")
    
    print("\nExample for t = -0.25:")
    print("The interval for a_i is [-1, -0.25]. The interval for sums is [-2, -0.5].")
    print("If a0 = -1 and a2 = -1, their sum is -2.")
    print("The required sum for a1+a3 is 1/(-2) = -0.5.")
    print("This sum can be formed by a1 = -0.25 and a3 = -0.25, which are in [-1, -0.25].")
    
    print("\nFinal Bounds:")
    print(f"{lower_bound} {upper_bound}")


find_bounds_of_t()