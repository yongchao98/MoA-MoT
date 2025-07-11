def find_bounds_of_t():
    """
    This function finds the lower and upper bounds of t by following a logical deduction.
    """
    print("This program finds the lower and upper bounds for a real number t based on a given condition.")
    print("The condition is: For any a_0, a_2 in [-1, t], there exist a_1, a_3 in [-1, t] such that (a_0 + a_2)(a_1 + a_3) = 1.\n")

    print("Step 1: The range of sums.")
    print("Let S = a_i + a_j. Since a_i, a_j are in [-1, t], the range of S is [-2, 2*t].")
    print("The condition implies that for any value x in [-2, 2*t], its reciprocal 1/x must also be in this range.\n")

    print("Step 2: The non-zero condition.")
    print("The term a_0 + a_2 cannot be zero, otherwise the equation 0 = 1 is impossible.")
    print("This means the range [-2, 2*t] cannot contain 0. Since -2 < 0, we must have 2*t < 0, which means t < 0.")
    print("The interval [-1, t] must also be valid, so -1 <= t. Thus, -1 <= t < 0.\n")

    print("Step 3: The set inclusion condition.")
    print("For t < 0, the range of sums [-2, 2*t] contains only negative numbers.")
    print("The set of reciprocals {1/x | x in [-2, 2*t]} is the interval [1/(2*t), -1/2].")
    print("This leads to the requirement: [1/(2*t), -1/2] must be a subset of [-2, 2*t].\n")

    print("Step 4: Solving the inequalities.")
    print("This subset condition yields two inequalities:")
    print("  1) 1/(2*t) >= -2")
    print("  2) -1/2 <= 2*t\n")

    print("Solving inequality 1) 1/(2*t) >= -2:")
    print("  Given t < 0, we multiply by 2*t and reverse the inequality sign: 1 <= -2 * (2*t)")
    num_1, num_neg_4 = 1, -4
    print(f"  {num_1} <= {num_neg_4}*t")
    print("  Divide by -4 and reverse the sign again: t <= 1/(-4)")
    upper_bound = num_1 / num_neg_4
    print(f"  Result 1: t <= {upper_bound}\n")

    print("Solving inequality 2) -1/2 <= 2*t:")
    print("  Divide by 2: t >= (-1/2) / 2")
    num_neg_half, num_2 = -0.5, 2
    lower_bound = num_neg_half / num_2
    print(f"  Result 2: t >= {lower_bound}\n")
    
    print("Step 5: Final conclusion.")
    print(f"Combining t <= {upper_bound} and t >= {lower_bound}, the only possible value for t is {lower_bound}.")
    print("Therefore, the lower and upper bounds of t are the same.\n")

    # Output the final answer
    print(f"{lower_bound} {upper_bound}")

find_bounds_of_t()