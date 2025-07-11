def solve_bounds():
    """
    This function solves for the lower and upper bounds of t by following a logical derivation.
    """
    print("Step-by-step derivation to find the bounds of t:")
    print("Let the interval for a_i be I = [-1, t].")
    print("Let S be the set of possible sums a_i + a_j. The range of S is [-2, 2t].")
    print("The given requirement is: For any s_0 in S, there must exist s_1 in S such that s_0 * s_1 = 1.")
    print("-" * 30)

    print("Step 1: Avoid the s_0 = 0 case.")
    print("The equation s_0 * s_1 = 1 has no solution if s_0 = 0.")
    print("The sum s_0 can be 0 if the range S = [-2, 2t] contains 0.")
    print("This occurs if 2t >= 0, which means t >= 0.")
    print("If t >= 0, we can choose a_0 = 0 and a_2 = 0. Their sum is 0, and the condition fails.")
    print("Therefore, we must have t < 0.")
    print("-" * 30)

    print("Step 2: Analyze the condition for t < 0.")
    print("If t < 0, the range S = [-2, 2t] contains only negative numbers.")
    print("The condition becomes: for any s in S, 1/s must also be in S.")
    print("This means the set {1/s | s in S} must be a subset of S.")
    print("The set S is the interval [-2, 2t].")
    print("The set of reciprocals {1/s | s in S} is the interval [1/(2t), 1/(-2)], which is [1/(2t), -0.5].")
    print("So, we require [1/(2t), -0.5] to be a subset of [-2, 2t].")
    print("-" * 30)

    print("Step 3: Formulate and solve the inequalities.")
    print("This requirement gives two inequalities:")
    print("  1. Lower bound inequality: -2 <= 1/(2t)")
    print("     Since t < 0, 2t is negative. Multiplying by 2t reverses the inequality sign.")
    print("     -2 * (2t) >= 1  =>  -4t >= 1  =>  t <= -1/4")
    print("\n  2. Upper bound inequality: -0.5 <= 2t")
    print("     Dividing by 2 gives: -0.25 <= t  =>  t >= -1/4")
    print("-" * 30)

    print("Step 4: Combine the results.")
    print("We have two conditions for t: t <= -1/4 and t >= -1/4.")
    print("The only value of t that satisfies both is t = -1/4.")
    
    lower_bound = -0.25
    upper_bound = -0.25
    
    print("\nFinal Conclusion:")
    print(f"The final equation for t is: t = {lower_bound}")
    print(f"This means the set of valid t is a single point.")
    print(f"The lower bound for t is {lower_bound}.")
    print(f"The upper bound for t is {upper_bound}.")
    
    print("\nAnswer:")
    print(f"{lower_bound} {upper_bound}")

solve_bounds()