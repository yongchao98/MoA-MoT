def solve_bounds():
    """
    This function solves for the lower and upper bounds of t based on the problem's requirements.
    The reasoning is explained through comments and print statements.
    """
    print("Step-by-step derivation of the bounds for t:")

    # The problem can be summarized as follows:
    # Let S = [-1, t] and R = {a+b | a, b in S} = [-2, 2t].
    # The condition is that for any x in R, 1/x must also be in R.

    # This condition leads to the set of reciprocals of R being a subset of R.
    # The set of reciprocals for x in [-2, 2t] (where t must be < 0 to avoid division by zero)
    # is [1/(2t), -1/2].
    print("The condition translates to the interval inclusion: [1/(2*t), -1/2] is a subset of [-2, 2*t].")

    # This inclusion gives two inequalities. Let's represent the numbers in these inequalities.
    # Inequality 1: -2 <= 1/(2t)
    # The numbers involved are -2, 1, 2.
    num_neg_2 = -2
    num_1 = 1
    num_2 = 2
    print(f"\n1. From the lower bounds: {num_neg_2} <= {num_1}/({num_2}*t)")
    print("   Assuming t < 0, we multiply by 2*t and flip the inequality sign.")
    # -2 * 2t >= 1 => -4t >= 1 => t <= -1/4
    t_upper_bound_val = -1/4
    print(f"   -4*t >= 1, which implies t <= {t_upper_bound_val}")
    
    # Inequality 2: -1/2 <= 2t
    # The numbers involved are -1, 2, 2.
    num_neg_1 = -1
    print(f"\n2. From the upper bounds: {num_neg_1}/{num_2} <= {num_2}*t")
    # t >= -1/4
    t_lower_bound_val = -1/4
    print(f"   Dividing by 2 gives t >= {t_lower_bound_val}")
    
    # The only value for t that satisfies both t <= -1/4 and t >= -1/4 is t = -1/4.
    solution_t = -0.25
    lower_bound = solution_t
    upper_bound = solution_t
    
    print(f"\nConclusion: The only value for t that satisfies both conditions is t = {solution_t}.")
    print("Therefore, the set of allowed values for t contains a single point.")
    
    # Print the final answer as requested
    print("\nThe lower and upper bounds for t are:")
    print(f"{lower_bound} {upper_bound}")

solve_bounds()