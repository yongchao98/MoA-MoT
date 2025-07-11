import math

def solve_random_walk_probability():
    """
    Calculates the probability for a 2D random walk to visit a set
    before leaving a disk.
    """
    # Parameters of the problem
    z0 = (0, 300)
    a1 = (0, 0)
    a2 = (2, 0)
    R = 1000
    r_eff = 1.0

    # Calculate required distances
    d = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)
    z0_a1_dist = math.sqrt((z0[0] - a1[0])**2 + (z0[1] - a1[1])**2)
    z0_a2_dist = math.sqrt((z0[0] - a2[0])**2 + (z0[1] - a2[1])**2)

    # Arguments for the logarithms in the formula
    # Numerator argument: R^2 / (|z0-a1| * |z0-a2|)
    num_arg = R**2 / (z0_a1_dist * z0_a2_dist)
    
    # Denominator argument: R^2 / (d * r_eff)
    den_arg = R**2 / (d * r_eff)

    # Calculate the numerator and the denominator
    numerator = math.log(num_arg)
    denominator = math.log(den_arg)

    # Calculate the final probability
    probability = numerator / denominator

    # Output the steps of the calculation as requested
    print("The probability P is calculated using the formula:")
    print("P = log(R^2 / (|z0-a1| * |z0-a2|)) / log(R^2 / (d * r_eff))")
    print("\nPlugging in the numbers:")
    print(f"R = {R}")
    print(f"|z0-a1| = {z0_a1_dist}")
    print(f"|z0-a2| = {z0_a2_dist:.4f}")
    print(f"d = {d}")
    print(f"r_eff = {r_eff}")
    print("\nEquation with values:")
    final_eq_str = (f"P = log({R**2:.0f} / ({z0_a1_dist:.0f} * {z0_a2_dist:.4f})) / "
                    f"log({R**2:.0f} / ({d:.0f} * {r_eff:.0f}))")
    print(final_eq_str)
    
    calc_str = (f"P = log({num_arg:.4f}) / log({den_arg:.0f}) = "
                f"{numerator:.4f} / {denominator:.4f} = {probability:.4f}")
    print(calc_str)

    # Final answer with three significant digits
    print(f"\nThe probability rounded to three significant digits is: {probability:.3g}")

solve_random_walk_probability()