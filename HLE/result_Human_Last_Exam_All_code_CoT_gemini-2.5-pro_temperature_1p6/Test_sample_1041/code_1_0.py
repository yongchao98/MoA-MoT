import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk visits a target set
    before leaving a large disk.
    """
    # Problem parameters
    start_pos = (0, 300)
    target1 = (0, 0)
    target2 = (2, 0)
    R = 1000.0

    # Calculate the distances needed for the formula
    dist_start_t1 = math.sqrt((start_pos[0] - target1[0])**2 + (start_pos[1] - target1[1])**2)
    dist_start_t2 = math.sqrt((start_pos[0] - target2[0])**2 + (start_pos[1] - target2[1])**2)
    dist_t1_t2 = math.sqrt((target1[0] - target2[0])**2 + (target1[1] - target2[1])**2)

    # Calculate the terms in the formula
    # u(z0) ~ (ln(R/|z0-t1|) + ln(R/|z0-t2|)) / ln(R^2/|t1-t2|)
    num1 = math.log(R / dist_start_t1)
    num2 = math.log(R / dist_start_t2)
    numerator = num1 + num2
    
    denominator = math.log(R**2 / dist_t1_t2)
    
    probability = numerator / denominator

    # Output the explanation and calculation steps
    print("The probability P is approximated by the formula:")
    print("P = [ln(R / |z_start - t1|) + ln(R / |z_start - t2|)] / [ln(R^2 / |t1 - t2|)]")
    print("\nWhere:")
    print(f"  R (disk radius) = {R}")
    print(f"  z_start (start position) = {start_pos}")
    print(f"  t1 (target 1) = {target1}")
    print(f"  t2 (target 2) = {target2}")
    print("\nCalculated distances:")
    print(f"  |z_start - t1| = {dist_start_t1:.4f}")
    print(f"  |z_start - t2| = {dist_start_t2:.4f}")
    print(f"  |t1 - t2| = {dist_t1_t2:.4f}")
    print("\nPlugging values into the formula:")
    
    final_eq = (f"P = (ln({R:.1f} / {dist_start_t1:.4f}) + ln({R:.1f} / {dist_start_t2:.4f}))"
                f" / ln({R**2:.1f} / {dist_t1_t2:.4f})")
    
    print(final_eq)
    
    final_calc = (f"P = ({num1:.4f} + {num2:.4f}) / {denominator:.4f}"
                  f" = {numerator:.4f} / {denominator:.4f}")

    print(final_calc)
    print(f"\nThe final probability is approximately: {probability:.3f}")


solve_random_walk_probability()