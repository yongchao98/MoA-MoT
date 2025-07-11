import math

def calculate_optimal_steps(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        kappa (float): The condition number of the function (M/mu).
                       It's assumed that mu=1 and M=kappa.

    Returns:
        tuple: A tuple containing the two optimal learning rates (gamma_1, gamma_2).
    """
    if kappa < 1:
        raise ValueError("Condition number kappa must be >= 1")
    if kappa == 1:
        return (1.0, 1.0)

    # The optimal polynomial gives the sum and product of the stepsizes
    # Denominator D = kappa^2 + 6*kappa + 1
    D = kappa**2 + 6 * kappa + 1
    
    # Sum of stepsizes: gamma_1 + gamma_2 = 8*(kappa+1) / D
    sum_gamma = 8 * (kappa + 1) / D
    
    # Product of stepsizes: gamma_1 * gamma_2 = 8 / D
    prod_gamma = 8 / D

    # The stepsizes are the roots of the quadratic equation:
    # x^2 - sum_gamma * x + prod_gamma = 0
    # We solve it using the quadratic formula: x = (-b +/- sqrt(b^2-4ac))/(2a)
    # Here, a=1, b=-sum_gamma, c=prod_gamma
    
    discriminant = sum_gamma**2 - 4 * prod_gamma
    
    gamma_1 = (sum_gamma - math.sqrt(discriminant)) / 2
    gamma_2 = (sum_gamma + math.sqrt(discriminant)) / 2
    
    return gamma_1, gamma_2

# --- Main execution ---
# Let's use kappa = 10 as an example value for the condition number M.
kappa = 10

# Calculate the optimal step sizes
gamma_1, gamma_2 = calculate_optimal_steps(kappa)

# Print the final result in the requested format
# The final equation is x^2 - (gamma_1+gamma_2)x + gamma_1*gamma_2 = 0
# and the solutions are gamma_1 and gamma_2
print("For a condition number kappa = {}".format(kappa))
print("The best choice for the pair (gamma_1, gamma_2) is:")
# We still need to output each number in the final equation!
# Let's derive the simplified closed-form expressions for printing
# D = kappa^2 + 6*kappa + 1 = 100 + 60 + 1 = 161
# S = 8*(10+1)/161 = 88/161
# P = 8/161
# The equation for the step-sizes is: x^2 - (88/161)x + 8/161 = 0
print("gamma_1 = (4 * ({} + 1) - 2 * sqrt(2) * ({} - 1)) / ({}**2 + 6 * {} + 1)".format(kappa, kappa, kappa, kappa))
print("        = (4 * ({}) - 2 * 1.4142 * ({})) / ({})".format(kappa + 1, kappa - 1, kappa**2 + 6 * kappa + 1))
print("        = {:.8f}".format(gamma_1))
print("")
print("gamma_2 = (4 * ({} + 1) + 2 * sqrt(2) * ({} - 1)) / ({}**2 + 6 * {} + 1)".format(kappa, kappa, kappa, kappa))
print("        = (4 * ({}) + 2 * 1.4142 * ({})) / ({})".format(kappa + 1, kappa - 1, kappa**2 + 6 * kappa + 1))
print("        = {:.8f}".format(gamma_2))
