import math

def calculate_alpha_from_psd():
    """
    Calculates the prompt neutron decay constant (alpha) from two
    Auto-Power Spectral Density (APSD) measurements.

    The model for the APSD is APSD(omega) = C / (alpha^2 + omega^2).
    Given two points (omega1, P1) and (omega2, P2), we can solve for alpha.
    The resulting formula is:
    alpha^2 = (P1*omega1^2 - P2*omega2^2) / (P2 - P1)
    """

    # --- Given values ---
    # Measurement 1:
    omega1 = 10.0  # rad/s
    P1 = 0.8     # APSD value at omega1 (arbitrary units)

    # Measurement 2:
    omega2 = 100.0 # rad/s
    P2 = 0.1      # APSD value at omega2 (arbitrary units)

    print("Problem: Calculate the prompt neutron decay constant, alpha, from two PSD measurements.")
    print("-" * 70)
    print(f"Given Measurement 1: Power (P1) = {P1}, Frequency (omega1) = {omega1} rad/s")
    print(f"Given Measurement 2: Power (P2) = {P2}, Frequency (omega2) = {omega2} rad/s")
    print("-" * 70)

    # Calculate the terms for the final equation
    numerator_term1 = P1 * omega1**2
    numerator_term2 = P2 * omega2**2
    denominator_term = P2 - P1

    # Check for division by zero, which happens if P1 == P2
    if denominator_term == 0:
        if numerator_term1 == numerator_term2:
            print("P1 equals P2 and the numerator is zero. Alpha cannot be determined from these points.")
        else:
            print("P1 equals P2, leading to division by zero. Cannot calculate alpha.")
        return

    # Calculate alpha squared
    alpha_squared = (numerator_term1 - numerator_term2) / denominator_term

    # Print the equation with the numbers plugged in, as requested
    print("The final equation for alpha^2 is: (P1 * omega1^2 - P2 * omega2^2) / (P2 - P1)")
    print("\nSubstituting the given values into the equation:")
    # The f-string formatting shows each number used in the calculation.
    print(f"alpha^2 = ({P1} * {omega1}^2 - {P2} * {omega2}^2) / ({P2} - {P1})")
    
    print("\nCalculating the terms:")
    print(f"alpha^2 = ({numerator_term1:.1f} - {numerator_term2:.1f}) / ({denominator_term:.1f})")
    
    numerator_result = numerator_term1 - numerator_term2
    print(f"alpha^2 = ({numerator_result:.1f}) / ({denominator_term:.1f})")

    if alpha_squared < 0:
        print("\nCalculation resulted in a negative value for alpha^2, which is physically unrealistic.")
        print(f"Result: alpha^2 = {alpha_squared:.4f}")
    else:
        # Calculate the final alpha
        alpha = math.sqrt(alpha_squared)
        print("\nFinal Result:")
        print(f"alpha^2 = {alpha_squared:.4f} s^-2")
        print(f"The prompt neutron decay constant, alpha = {alpha:.4f} s^-1")


# Execute the function
if __name__ == "__main__":
    calculate_alpha_from_psd()