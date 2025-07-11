def find_bounds_and_print():
    """
    This script finds the lower and upper bounds for the real number t
    by codifying the logical derivation explained above.
    """
    # Plan:
    # 1. The problem requires the range of sums K=[-2, 2t] to be closed under the reciprocal operation.
    # 2. This implies that for K, the set of reciprocals f(K) must be a subset of K.
    # 3. This requirement translates to two inequalities for t:
    #    - t <= -1/4 (sets the upper bound)
    #    - t >= -1/4 (sets the lower bound)
    # 4. We solve the boundary equations to find the exact values for the bounds.

    print("To find the lower bound, we solve the boundary condition from the inequality -1/2 <= 2*t.")
    # The boundary is at -1/2 = 2*t.
    num1_l, den1_l = -1, 2
    num2_l = 2
    lower_bound = (num1_l / den1_l) / num2_l
    print(f"The final equation for the lower bound is: ({num1_l}/{den1_l}) / {num2_l} = {lower_bound}")
    
    print("\nTo find the upper bound, we solve the boundary condition from the inequality -2 <= 1/(2*t).")
    # The boundary is at -2 = 1/(2*t), which simplifies to -4*t = 1, so t = 1/(-4).
    num1_u = 1
    num2_u = -4
    upper_bound = num1_u / num2_u
    print(f"The final equation for the upper bound is: {num1_u} / ({num2_u}) = {upper_bound}")

    print(f"\n{lower_bound} {upper_bound}")

if __name__ == "__main__":
    find_bounds_and_print()