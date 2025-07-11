import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem near the ground.

    The problem is solved using the point vortex model and the mirror image method.
    The final lift ratio L1/L2 is given by (1+K)/(1-K), where K is an
    interaction factor dependent on the geometry.
    """
    # We can set the chord 'c' to any value since it will cancel out. Let's use 2.0 for simplicity.
    c = 2.0

    # Given geometric parameters
    s = 0.5 * c  # Separation s = 1/2 * c
    h = 0.5 * c  # Ride height h = c/2

    print("--- Problem Setup ---")
    print(f"The aerofoils are modeled as point vortices in a freestream.")
    print(f"Ground effect is handled by the mirror image method.")
    print(f"Chord length, c = {c:.1f}")
    print(f"Horizontal separation, s = c/2 = {s:.1f}")
    print(f"Ride height, h = c/2 = {h:.1f}")
    print("\n--- Calculation Steps ---")

    # Calculate the interaction factor K
    # K = c/(2*s) - c*s/(2*(s**2 + 4*h**2))
    term1 = c / (2 * s)
    term2_numerator = c * s
    term2_denominator = 2 * (s**2 + 4 * h**2)
    term2 = term2_numerator / term2_denominator
    K = term1 - term2

    print("1. Calculate the interaction factor, K:")
    print(f"K = c/(2*s) - c*s/(2*(s^2 + 4*h^2))")
    print(f"K = {c:.1f}/(2*{s:.1f}) - ({c:.1f}*{s:.1f})/(2*({s:.1f}^2 + 4*{h:.1f}^2))")
    print(f"K = {term1:.2f} - {term2_numerator:.2f}/{term2_denominator:.2f}")
    print(f"K = {term1:.2f} - {term2:.2f} = {K:.2f}")
    print("")

    # Calculate the lift ratio L1/L2
    # L1/L2 = (1 + K) / (1 - K)
    lift_ratio = (1 + K) / (1 - K)

    print("2. Calculate the lift ratio, L1/L2:")
    print(f"L1/L2 = (1 + K) / (1 - K)")
    print(f"L1/L2 = (1 + {K:.2f}) / (1 - {K:.2f})")
    print(f"L1/L2 = {1 + K:.2f} / {1 - K:.2f}")
    print(f"L1/L2 = {lift_ratio:.2f}")
    print("\n--- Final Answer ---")
    print(f"The calculated lift ratio L1/L2 is: {lift_ratio}")

if __name__ == "__main__":
    calculate_lift_ratio()