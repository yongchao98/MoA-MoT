def solve_bounds():
    """
    This function follows the step-by-step plan to determine the bounds of t.
    """
    # Step 1 & 2: Define variables and ranges
    # a_i are in [-1, t]
    # S = a_i + a_j is in [-2, 2t]
    # Condition: S1 * S2 = 1

    # Step 3, 4, 5: Handle the zero-sum case and constrain t
    # To avoid S1 = 0, the interval [-2, 2t] must not contain 0.
    # Since -2 is a fixed lower bound, the interval must be entirely negative.
    # This requires 2t < 0, which means t < 0.
    # For the interval [-1, t] to be valid, t >= -1.
    # So, we have established that -1 <= t < 0.
    print("Step 1: To prevent the sum a_0 + a_2 from being zero, we must have t < 0.")

    # Step 6 & 7: Set up and solve the main inequality
    # The condition is that for any s in [-2, 2t], 1/s must also be in [-2, 2t].
    # The image of [-2, 2t] under f(s)=1/s is [1/(2t), -1/2].
    # So, we need [1/(2t), -1/2] to be a subset of [-2, 2t].
    # This gives two inequalities:
    # 1) -1/2 <= 2t  =>  t >= -1/4
    # 2) 1/(2t) >= -2 => 1 <= -4t (since t<0, 2t is negative) => t <= -1/4
    print("Step 2: The problem reduces to solving the inclusion [1/(2t), -1/2] ⊆ [-2, 2t].")
    print("This yields two inequalities for t:")
    print("  a) From the upper bounds: -1/2 <= 2t  =>  t >= -1/4")
    print("  b) From the lower bounds: 1/(2t) >= -2  =>  t <= -1/4")

    # Step 8: Determine the final bounds
    # The conditions t >= -1/4 and t <= -1/4 mean t must be exactly -1/4.
    lower_bound = -0.25
    upper_bound = -0.25
    print("\nStep 3: The only value of t satisfying both inequalities is t = -0.25.")
    print(f"Therefore, the lower bound is {lower_bound} and the upper bound is {upper_bound}.")
    
    print("\nFinal Answer:")
    print(f"{lower_bound} {upper_bound}")

solve_bounds()