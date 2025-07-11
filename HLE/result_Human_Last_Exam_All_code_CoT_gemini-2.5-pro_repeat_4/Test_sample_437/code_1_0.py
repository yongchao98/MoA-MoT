import math

def evaluate_1s_integral(zeta):
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-type orbital.

    This function follows the analytical derivation step-by-step and prints
    the result for a given orbital exponent zeta. This integral is known
    as a one-electron integral and represents the electron-nucleus attraction energy
    (when multiplied by -Z, where Z is the nuclear charge).
    """
    print("### Evaluation of the integral <phi_i | 1/r | phi_j> for i=j=1s STO ###\n")

    # Step 1: Define the orbital and constants
    print("1. The normalized 1s Slater-Type Orbital (STO) is given by:")
    print(f"   phi_1s(r) = (zeta^3 / pi)^(1/2) * exp(-zeta*r)")
    print(f"   For this calculation, we will use a realistic exponent for Hydrogen, zeta = {zeta}\n")

    # Step 2: Set up the integral expression
    print("2. The integral to evaluate is I = Integral( phi_1s(r) * (1/r) * phi_1s(r) * dV ) over all space.")
    print("   Substituting the STO expression and dV = r^2 * sin(theta) * dr * d(theta) * d(phi):")
    print("   I = Integral( (zeta^3/pi) * exp(-2*zeta*r) * (1/r) * r^2 * sin(theta) * dr * d(theta) * d(phi) )\n")

    # Step 3: Separate into angular and radial parts
    print("3. The integral is separable into three parts: Normalization, Angular, and Radial.")
    print("   I = (zeta^3 / pi) * [Angular Integral] * [Radial Integral]\n")
    norm_squared = zeta**3 / math.pi

    # Step 4: Evaluate the angular part
    print("4. Evaluating the angular part:")
    print("   Angular Integral = Integral from 0 to 2pi d(phi) * Integral from 0 to pi sin(theta) d(theta)")
    angular_part = 4 * math.pi
    print(f"   The result of the angular integration is exactly 4*pi.\n")

    # Step 5: Evaluate the radial part
    print("5. Evaluating the radial part:")
    print("   The integral simplifies to: Integral from 0 to infinity of (r * exp(-2*zeta*r)) dr")
    print("   We use the standard integral formula: Integral(x^n * exp(-a*x) dx) = n! / a^(n+1)")
    print(f"   For our integral, n=1 and a = 2*zeta = 2 * {zeta} = {2*zeta:.2f}")
    radial_part = 1 / ((2 * zeta)**2)
    print(f"   The result of the radial integration is 1 / (2*zeta)^2.\n")

    # Step 6: Combine all parts and show the final equation with numbers
    print("6. Combining all parts to get the final result:")
    print("   I = (Normalization_Constant^2) * (Angular Part) * (Radial Part)")
    
    final_result = norm_squared * angular_part * radial_part

    print("\n   Final Equation (symbolic):")
    print("   I = (zeta^3 / pi) * (4*pi) * (1 / (4*zeta^2))")
    print("   This expression simplifies analytically to I = zeta.")

    print("\n   Final Equation (with numerical values for each part):")
    print(f"   I = ({norm_squared:.4f}) * ({angular_part:.4f}) * ({radial_part:.4f})")
    
    print(f"\n   Multiplying these values: {norm_squared:.4f} * {angular_part:.4f} * {radial_part:.4f} = {final_result:.2f}")

    print(f"\n   As shown by both the analytical simplification and the numerical calculation,")
    print(f"   the final value of the integral is equal to zeta.")


# --- Main execution ---
# Use a realistic zeta value for the Hydrogen 1s orbital in an STO-3G basis set.
zeta_value = 1.24
evaluate_1s_integral(zeta_value)

print(f"\nFinal Answer: {zeta_value}")