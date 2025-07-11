import sympy

def find_t_bounds():
    """
    This function finds the lower and upper bounds of the real number t
    by following a step-by-step logical derivation from the problem statement.
    """
    # Let t be the real number we are looking for.
    # The variables a_0, a_1, a_2, a_3 are in the interval [-1, t].
    # Let S1 = a_0 + a_2 and S2 = a_1 + a_3.
    # The range of possible values for S1 and S2 is the interval I = [-1 + (-1), t + t] = [-2, 2t].

    # Step 1: Establish the necessary condition that t < 0.
    # The problem states (a_0 + a_2)(a_1 + a_3) = 1.
    # This implies that the sum a_0 + a_2 can never be zero. If it were, the equation would become 0 = 1,
    # which is impossible to satisfy by any choice of a_1 and a_3.
    # The range of a_0 + a_2 is the interval I = [-2, 2t].
    # For 0 not to be in I, the interval must be entirely on one side of 0.
    # Since the lower bound is -2, the interval must be entirely negative.
    # This requires the upper bound to be negative: 2t < 0, which implies t < 0.
    print("Step 1: Deduce the sign of t.")
    print("Let S = a_0 + a_2. The range of S is I = [-2, 2t].")
    print("If S could be 0, the equation (a_0+a_2)(a_1+a_3)=1 would be 0=1, which is impossible.")
    print("Therefore, 0 must not be in the range of S. This means 2t < 0, so t < 0.\n")

    # Step 2: Set up inequalities for t based on the core requirement.
    # The condition is that for any value s1 in I, there must exist s2 in I such that s1 * s2 = 1.
    # This means that for any s1 in I, its reciprocal 1/s1 must also be in I.
    # So, the set {1/x | x in I} must be a subset of I.
    # Since t < 0, the interval I = [-2, 2t] contains only negative numbers.
    # The function f(x) = 1/x is decreasing on this interval.
    # The image of I under f(x) is [1/(2t), 1/(-2)] = [1/(2t), -1/2].
    # For this image to be a subset of I, we need:
    #   1) The lower bound of the image >= the lower bound of I.
    #   2) The upper bound of the image <= the upper bound of I.
    print("Step 2: Set up inequalities for t.")
    print("The condition implies that for any s in I = [-2, 2t], 1/s must also be in I.")
    print("The set of reciprocals {1/s | s in [-2, 2t]} is the interval [1/(2t), -1/2].")
    print("So, we require [1/(2t), -1/2] to be a subset of [-2, 2t].")
    print("This gives two inequalities:")
    print("  1) -2 <= 1/(2t)")
    print("  2) -1/2 <= 2t\n")

    # Step 3: Solve the inequalities.
    # We use symbolic representation for clarity.
    t_sym = sympy.Symbol('t')
    
    # Solving inequality 1: -2 <= 1/(2t)
    # Since t < 0, 2t is negative. Multiplying by 2t reverses the inequality sign.
    # -2 * (2t) >= 1  =>  -4t >= 1  =>  t <= -1/4
    print("Step 3: Solve the inequalities.")
    print("Solving -2 <= 1/(2t) for t < 0 gives: t <= -1/4.")
    
    # Solving inequality 2: -1/2 <= 2t
    # Dividing by 2 (a positive number) does not change the inequality.
    # -1/4 <= t
    print("Solving -1/2 <= 2t gives: t >= -1/4.")
    print("Combining t <= -1/4 and t >= -1/4, the only possible value is t = -1/4.\n")

    # Step 4: State the final equation and bounds.
    # The only value for t that satisfies the conditions is t = -1/4.
    t_solution = -sympy.S(1)/4
    lower_bound = float(t_solution)
    upper_bound = float(t_solution)

    print("Step 4: State the final equation and bounds.")
    print(f"The final equation for t is: t = {t_solution}")
    print(f"The lower bound for t is {lower_bound}.")
    print(f"The upper bound for t is {upper_bound}.")
    
    # Final answer as requested by the user prompt.
    print("\nFinal Answer:")
    print(f"{lower_bound} {upper_bound}")

find_t_bounds()