import math

def calculate_frequency_correction(gamma):
    """
    Calculates the second-order nonlinear frequency correction (omega_2) for the
    dimensionless Rayleigh-Plesset equation.

    The formula derived using the Poincaré-Lindstedt method is:
    omega_2 = - (sqrt(3*gamma) * (24*gamma^2 + 15*gamma + 2)) / 16
    """

    # Print the explanation of the task
    print("This script calculates the coefficient of the first non-zero term in the nonlinear")
    print("correction to the bubble oscillation frequency, denoted as omega_2.")
    print(f"We will use a standard value for the polytropic index of air, gamma = {gamma}\n")

    # Calculate each part of the equation
    numerator_term1 = math.sqrt(3 * gamma)
    polynomial_term = 24 * gamma**2 + 15 * gamma + 2
    denominator = 16

    # Calculate the final result for omega_2
    omega_2 = - (numerator_term1 * polynomial_term) / denominator

    # Output the equation with the calculated numbers as requested
    print("The formula for the frequency correction omega_2 is:")
    print(f"omega_2 = - (sqrt(3 * {gamma}) * (24 * {gamma}^2 + 15 * {gamma} + 2)) / 16\n")
    print("Substituting the value of gamma:")
    print(f"Numerator Part 1 (sqrt(3 * gamma)): sqrt({3*gamma:.2f}) = {numerator_term1:.4f}")
    print(f"Numerator Part 2 (Polynomial): 24*{gamma**2:.2f} + 15*{gamma:.1f} + 2 = {24*gamma**2:.2f} + {15*gamma:.1f} + 2 = {polynomial_term:.4f}")
    print(f"Denominator: {denominator}\n")

    print("The final equation is:")
    print(f"omega_2 = - ({numerator_term1:.4f} * {polynomial_term:.4f}) / {denominator}")
    print(f"omega_2 = - ({numerator_term1 * polynomial_term:.4f}) / {denominator}\n")

    print("The numerical value for the coefficient omega_2 is:")
    print(f"{omega_2:.4f}")
    
    return omega_2

if __name__ == '__main__':
    # Assume a value for the polytropic index gamma, e.g., for diatomic gas (air).
    gamma_value = 1.4
    final_answer = calculate_frequency_correction(gamma_value)
    # The final answer in the required format is the numerical value of omega_2
    print(f"\n<<<{final_answer:.4f}>>>")
