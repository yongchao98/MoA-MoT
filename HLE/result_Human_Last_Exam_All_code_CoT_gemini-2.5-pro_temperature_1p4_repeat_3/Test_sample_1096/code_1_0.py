import math

def calculate_optimal_gammas(kappa):
    """
    Calculates and prints the optimal step sizes (gamma_1, gamma_2)
    for two-step gradient descent based on the condition number kappa.

    Args:
        kappa (float): The condition number M/mu, must be > 1.
    """
    if kappa <= 1:
        print("Error: Kappa must be greater than 1.")
        return

    # The formulas are derived from Chebyshev polynomials.
    # The optimal gamma values are the roots of a quadratic equation
    # x^2 - sum_gamma * x + prod_gamma = 0
    # where sum_gamma and prod_gamma depend on kappa.

    print(f"Calculating optimal step sizes for kappa = {kappa}:")

    # Common denominator in the expressions for gammas
    denominator = kappa**2 + 6*kappa + 1
    
    # Sum and product of the optimal gammas
    sum_gamma = 8 * (kappa + 1) / denominator
    prod_gamma = 8 / denominator
    
    print("\nThe optimal step sizes (gamma_1, gamma_2) are the roots of the quadratic equation:")
    print(f"x^2 - ({sum_gamma:.6f})*x + ({prod_gamma:.6f}) = 0")


    # We can solve this quadratic equation directly. The roots are:
    # gamma = (sum_gamma +/- sqrt(sum_gamma^2 - 4*prod_gamma)) / 2
    # which simplifies to the expressions below.

    # Coefficients used in the direct formulas for gamma_1 and gamma_2
    c1 = 4 + 2 * math.sqrt(2)
    c2 = 4 - 2 * math.sqrt(2)

    # Numerators for gamma_1 and gamma_2
    numerator_1 = c1 * kappa + c2
    numerator_2 = c2 * kappa + c1
    
    # Calculate final gamma values
    gamma_1 = numerator_1 / denominator
    gamma_2 = numerator_2 / denominator

    # Output the detailed calculation for transparency
    # This follows the request "output each number in the final equation"
    print("\n--- Detailed Calculation ---")
    print("\nThe general formulas for the step sizes are:")
    print("gamma_1 = ((4 + 2*sqrt(2))*kappa + (4 - 2*sqrt(2))) / (kappa^2 + 6*kappa + 1)")
    print("gamma_2 = ((4 - 2*sqrt(2))*kappa + (4 + 2*sqrt(2))) / (kappa^2 + 6*kappa + 1)")
    
    print("\nFor gamma_1:")
    print(f"  Numerator   = ({c1:.6f} * {kappa}) + {c2:.6f} = {numerator_1:.6f}")
    print(f"  Denominator = {kappa}^2 + 6*{kappa} + 1 = {denominator:.6f}")
    print(f"  gamma_1 = {numerator_1:.6f} / {denominator:.6f} = {gamma_1:.8f}")

    print("\nFor gamma_2:")
    print(f"  Numerator   = ({c2:.6f} * {kappa}) + {c1:.6f} = {numerator_2:.6f}")
    print(f"  Denominator = {kappa}^2 + 6*{kappa} + 1 = {denominator:.6f}")
    print(f"  gamma_2 = {numerator_2:.6f} / {denominator:.6f} = {gamma_2:.8f}")
    
    print("\nNote: The order of gamma_1 and gamma_2 is interchangeable.")


# --- User-configurable Section ---
# Set the condition number kappa (M/mu). A common example value is 10.
# You can change this value to suit your specific function.
kappa_value = 10.0
# --- End of User-configurable Section ---

calculate_optimal_gammas(kappa_value)