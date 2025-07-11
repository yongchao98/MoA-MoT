def solve_bounds():
    """
    This function solves for the lower and upper bounds of t based on the problem's requirement.
    The reasoning is implemented step-by-step.
    """

    # The problem requires finding the bounds for a real number t.
    # The condition is: For any a_0, a_2 in [-1, t], there exist a_1, a_3 in [-1, t]
    # such that (a_0 + a_2)(a_1 + a_3) = 1.

    # Let R be the range of the sum a_i + a_j, where a_i, a_j are in [-1, t].
    # The minimum sum is -1 + (-1) = -2.
    # The maximum sum is t + t = 2t.
    # So, R = [-2, 2t].

    # The condition can be rephrased: For any x in R, 1/x must also be in R.
    # (Excluding x=0, which we handle next).

    # Step 1: The sum cannot be zero.
    # If x = a_0 + a_2 = 0, the equation 0 * y = 1 has no solution.
    # The range R = [-2, 2t] contains 0 if 2t >= 0, i.e., t >= 0.
    # To prevent the sum from ever being 0, the range R must not contain 0.
    # This implies 2t < 0, so t < 0.
    # Also, for the interval [-1, t] to be valid, t >= -1.
    # So, t must be in the range [-1, 0).

    # Step 2: For t in [-1, 0), find the condition for R.
    # R = [-2, 2t]. This is an interval of negative numbers.
    # The set of reciprocals is R' = {1/x for x in R}.
    # Since f(x)=1/x is decreasing for x<0, the reciprocal set is [1/(2t), 1/(-2)].
    # R' = [1/(2t), -1/2].

    # The condition is R' must be a subset of R.
    # [1/(2t), -1/2] must be a subset of [-2, 2t].

    # This leads to two inequalities:
    # 1) 1/(2t) >= -2  => 1 <= -4t (since t<0, flip inequality) => t <= -1/4
    # 2) -1/2 <= 2t   => t >= -1/4

    # For both inequalities to hold, t must be exactly -1/4.
    t_solution = -1.0 / 4.0

    # The set of valid t is {-1/4}.
    lower_bound = t_solution
    upper_bound = t_solution

    # Print the final answer
    print(f"{lower_bound} {upper_bound}")

solve_bounds()