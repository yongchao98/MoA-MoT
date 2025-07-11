import math

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem in ground effect.

    This function follows these steps:
    1.  Sets up the geometry based on the problem statement.
    2.  Calculates the interaction coefficient 'K' which encapsulates the geometric
        effects of tandem separation and ground proximity.
    3.  Solves for the lift ratio using the derived formula L1/L2 = (1 + K) / (1 - K).
    4.  Prints the steps of the calculation and the final result.
    """

    # Step 1: Define the geometry.
    # We can set the chord 'c' to 1.0, as it will cancel out in the final ratio.
    c = 1.0
    s = 0.5 * c  # Horizontal separation
    h = 0.5 * c  # Ride height

    print("--- Aerodynamic Calculation for Tandem Aerofoils in Ground Effect ---")
    print("\n1. System Geometry:")
    print(f"   - Aerofoil Chord, c = {c:.1f}")
    print(f"   - Separation, s = 1/2 * c = {s:.2f}")
    print(f"   - Ride Height, h = 1/2 * c = {h:.2f}")

    # Step 2: Calculate the interaction coefficient K.
    # The interaction coefficient K is derived from potential flow theory and is given by:
    # K = (c/2) * [1/s - s / (s^2 + 4*h^2)]
    print("\n2. Calculating the Interaction Coefficient (K):")
    print("   The formula for K is: K = (c/2) * [1/s - s / (s^2 + 4*h^2)]")

    # Breaking down the calculation for clarity
    term1_inv_s = 1 / s
    term2_numerator = s
    term2_denominator = s**2 + 4 * h**2
    term2_fraction = term2_numerator / term2_denominator

    print(f"   Substituting values:")
    print(f"   1/s = 1 / {s:.2f} = {term1_inv_s:.4f}")
    print(f"   s / (s^2 + 4*h^2) = {s:.2f} / ({s:.2f}^2 + 4*{h:.2f}^2) = {term2_numerator:.2f} / {term2_denominator:.2f} = {term2_fraction:.4f}")

    K = (c / 2) * (term1_inv_s - term2_fraction)
    print(f"\n   K = ({c:.1f}/2) * [{term1_inv_s:.4f} - {term2_fraction:.4f}]")
    print(f"   K = {c/2:.1f} * [{term1_inv_s - term2_fraction:.4f}] = {K:.4f}")


    # Step 3: Calculate the final lift ratio L1/L2.
    # From the system of equations for the vortex circulations, we get:
    # L1/L2 = Γ1/Γ2 = (1 + K) / (1 - K)
    print("\n3. Calculating the Lift Ratio (L1/L2):")
    print("   The formula for the lift ratio is: L1/L2 = (1 + K) / (1 - K)")

    lift_ratio = (1 + K) / (1 - K)
    numerator = 1 + K
    denominator = 1 - K

    print(f"   Substituting K = {K:.2f}:")
    # This part prints the final equation with numbers, as requested
    print(f"   L1 / L2 = (1 + {K:.2f}) / (1 - {K:.2f})")
    print(f"           = {numerator:.2f} / {denominator:.2f}")
    print(f"           = {lift_ratio:.1f}")

    print("\n--- End of Calculation ---")
    return lift_ratio

# Execute the function
final_ratio = calculate_lift_ratio()
# The final answer is wrapped in <<<>>> as per instructions.
# print(f"\nFinal Answer: <<< {final_ratio:.1f} >>>") # This line is for internal check

if __name__ == '__main__':
    pass # The function call above runs the script.