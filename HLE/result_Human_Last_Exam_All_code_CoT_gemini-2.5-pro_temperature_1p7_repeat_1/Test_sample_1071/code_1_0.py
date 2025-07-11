def find_t_bounds():
    """
    This script finds the lower and upper bounds for the real number t
    based on the given mathematical condition.
    """
    print("Step 1: Define the sums and their ranges.")
    print("Let S_sum = a_i + a_j, where a_i and a_j are in the interval [-1, t].")
    print("The smallest possible sum is (-1) + (-1) = -2.")
    print("The largest possible sum is t + t = 2t.")
    print("So, the range of possible sums, let's call it R, is [-2, 2t].")
    print("The condition (a_0 + a_2)(a_1 + a_3) = 1 means that for any value x in R,")
    print("the value 1/x must also be in R, since a_1, a_3 can also be chosen from [-1, t].")
    print("-" * 30)

    print("Step 2: Handle the possibility of division by zero.")
    print("The expression 1/x is undefined if x = 0.")
    print("If 0 is in the range R = [-2, 2t], we could pick a_0 and a_2 such that a_0 + a_2 = 0.")
    print("In this case, the equation 0 * (a_1 + a_3) = 1 would have no solution, and the requirement would fail.")
    print("To prevent this, the range R = [-2, 2t] must not contain 0.")
    print("This implies either the entire range is positive (-2 > 0, which is false) or the entire range is negative (2t < 0).")
    print("Therefore, we must have 2t < 0, which means t < 0.")
    print("-" * 30)

    print("Step 3: Solve the main inequality for t < 0.")
    print("With t < 0, the range R = [-2, 2t] consists of only negative numbers.")
    print("Our condition is: for any x in [-2, 2t], 1/x must also be in [-2, 2t].")
    print("Let's find the range of 1/x. Since x is negative, 1/x is also negative.")
    print("The function f(x) = 1/x is increasing for x < 0.")
    print("So, the range of 1/x for x in [-2, 2t] is [1/(-2), 1/(2t)] = [-0.5, 1/(2t)].")
    print("This new range must be a subset of the original range R, so [-0.5, 1/(2t)] must be a subset of [-2, 2t].")
    print("This gives us two conditions:")
    print("  1) -2 <= -0.5  (This is always true).")
    print("  2) 1/(2t) <= 2t")
    print("-" * 30)
    
    print("Step 4: Solve the final inequality.")
    print("We solve: 1/(2*t) <= 2*t")
    print("Since t < 0, 2*t is a negative number. When we multiply both sides by 2*t, we must flip the inequality sign.")
    print("This gives: 1 >= (2*t)*(2*t)")
    print("Equation: 1 >= 4*t^2")
    print("Dividing by 4: t^2 <= 1/4")
    print("Taking the square root: |t| <= sqrt(1/4), so |t| <= 0.5")
    print("Because we know t is negative (t < 0), |t| is equal to -t.")
    print("So, -t <= 0.5, which means t >= -0.5.")
    print("-" * 30)

    print("Step 5: Combine all constraints on t.")
    print("The initial problem requires the interval [-1, t] to be valid, so t >= -1.")
    print("From Step 2, we must have t < 0.")
    print("From Step 4, we must have t >= -0.5.")
    print("Combining all three (t >= -1, t < 0, and t >= -0.5), we get the final range for t: -0.5 <= t < 0.")
    print("-" * 30)
    
    lower_bound = -0.5
    upper_bound = 0
    
    print("Conclusion:")
    print(f"The set of all valid values for t is the interval [{lower_bound}, {upper_bound}).")
    print("The lower bound is the minimum value in this set.")
    print("The upper bound is the supremum (least upper bound) of this set.")
    print("\nFinal Answer:")
    print(f"{lower_bound} {upper_bound}")

find_t_bounds()