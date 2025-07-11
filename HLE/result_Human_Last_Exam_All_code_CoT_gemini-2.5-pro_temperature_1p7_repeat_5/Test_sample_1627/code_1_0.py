import sys

def calculate_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # For calculation purposes, let's assume the inner radius 'a' is 1.0.
    # The actual value doesn't matter since we are calculating a ratio.
    a = 1.0

    # The problem states the outer radius 'b' is twice the inner radius 'a'.
    b = 2.0 * a

    print("In elasticity theory, the maximum tangential stress (σ_t_max) in a thick-walled cylinder")
    print("under internal pressure (P_i) is given by the formula:")
    print("σ_t_max = P_i * (b² + a²) / (b² - a²)")
    print("where 'a' is the inner radius and 'b' is the outer radius.\n")

    print("We need to find the ratio: σ_t_max / P_i")
    print("Ratio = (b² + a²) / (b² - a²)\n")

    # Given condition: b = 2a. Let's use a=1 and b=2.
    a_squared = a**2
    b_squared = b**2

    print("Substituting the values (a=1.0, b=2.0) into the ratio equation:")
    # Print the equation with numbers
    # We use f-strings with specified precision for clear output.
    print(f"Ratio = ({b:.1f}² + {a:.1f}²) / ({b:.1f}² - {a:.1f}²)")
    print(f"Ratio = ({b_squared:.1f} + {a_squared:.1f}) / ({b_squared:.1f} - {a_squared:.1f})")

    # Perform the calculation
    numerator = b_squared + a_squared
    denominator = b_squared - a_squared
    ratio = numerator / denominator

    print(f"Ratio = {numerator:.1f} / {denominator:.1f}")
    
    # In python 3, division of integers produces float.
    # To represent as a fraction, we can show numerator and denominator
    # For example using `fractions` module, or simply print it.
    print(f"\nFinal Ratio = 5 / 3")
    print(f"As a decimal, the ratio is approximately: {ratio:.4f}")

if __name__ == "__main__":
    calculate_stress_ratio()
