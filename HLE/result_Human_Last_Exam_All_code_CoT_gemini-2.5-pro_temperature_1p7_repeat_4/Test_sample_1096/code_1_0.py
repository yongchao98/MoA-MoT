import math

def calculate_optimal_step_sizes(kappa):
    """
    Calculates the optimal step sizes (gamma_1, gamma_2) for a two-step gradient descent method.

    Args:
        kappa (float): The condition number of the function (M/m). Must be > 1.

    Returns:
        tuple: A tuple containing the two optimal step sizes (gamma_1, gamma_2).
    """
    if kappa <= 1:
        raise ValueError("Condition number kappa must be greater than 1.")

    # The reciprocals of the optimal step sizes are the roots of the optimal polynomial
    # derived from Chebyshev polynomials.
    # 1/gamma = (kappa+1)/2 +/- (kappa-1)/(2*sqrt(2))
    
    # Let's define the two roots
    root1 = (kappa + 1) / 2 + (kappa - 1) / (2 * math.sqrt(2))
    root2 = (kappa + 1) / 2 - (kappa - 1) / (2 * math.sqrt(2))
    
    # The step sizes are the reciprocals of these roots.
    # To maintain a consistent order, we can assign the smaller step size to gamma_1.
    # A larger root corresponds to a smaller step size.
    gamma_1 = 1 / root1
    gamma_2 = 1 / root2
    
    return gamma_1, gamma_2

# --- Example Usage ---
# You can change the value of kappa here.
# It is assumed that mu=1, M=kappa.
kappa = 10

# Calculate the optimal step sizes
gamma_1, gamma_2 = calculate_optimal_step_sizes(kappa)

# Print the final result in a formatted equation
print("The best choice for the pair (gamma_1, gamma_2) is:")
# We want to print the full equation with the numbers plugged in.
# gamma_1 formula
term1_g1 = (kappa+1)/2
term2_g1 = (kappa-1)/(2*math.sqrt(2))
print(f"gamma_1 = 1 / (({kappa}+1)/2 + ({kappa}-1)/(2*sqrt(2))) = 1 / ({term1_g1:.4f} + {term2_g1:.4f}) = {gamma_1:.4f}")

# gamma_2 formula
term1_g2 = (kappa+1)/2
term2_g2 = (kappa-1)/(2*math.sqrt(2))
print(f"gamma_2 = 1 / (({kappa}+1)/2 - ({kappa}-1)/(2*sqrt(2))) = 1 / ({term1_g2:.4f} - {term2_g2:.4f}) = {gamma_2:.4f}")
