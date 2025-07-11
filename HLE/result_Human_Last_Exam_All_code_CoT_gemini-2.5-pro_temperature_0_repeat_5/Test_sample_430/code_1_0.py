import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site based on a
    simplified multi-layer mean-field model.
    """
    # Dimensionless parameters given in the problem
    # mu_prime = mu / (k_B * T)
    mu_prime = 0.15
    # eps1_prime = epsilon_1 / (k_B * T)
    eps1_prime = 0.1

    # In the limit of infinite layers (k -> infinity), epsilon_l -> 0 and epsilon_2 -> 0.
    # The model simplifies to layer-by-layer Langmuir adsorption.

    # Calculate the conditional coverage of the first layer (theta_1)
    # theta_1 = 1 / (1 + exp(-(mu' + eps1')))
    theta_1 = 1 / (1 + math.exp(-(mu_prime + eps1_prime)))

    # Calculate the conditional coverage of subsequent layers (theta_2)
    # Since eps2' -> 0, the equation is theta_2 = 1 / (1 + exp(-mu'))
    theta_2 = 1 / (1 + math.exp(-mu_prime))

    # The average number of layers <k> is the sum of a geometric series:
    # <k> = theta_1 / (1 - theta_2)
    avg_layers = theta_1 / (1 - theta_2)

    # Print the final equation with numerical values
    print("The equation for the average number of layers <k> is:")
    print("<k> = θ₁ / (1 - θ₂)")
    print("   = [1 / (1 + exp(-(μ/kBT + ε₁/kBT)))] / [1 - 1 / (1 + exp(-μ/kBT))]")
    print(f"   = [1 / (1 + exp(-({mu_prime} + {eps1_prime})))] / [1 - 1 / (1 + exp(-{mu_prime}))]")
    print(f"   = {theta_1:.4f} / (1 - {theta_2:.4f})")
    
    # Print the final result
    print("\nThe calculated average number of adsorbed layers per site is:")
    print(f"{avg_layers}")
    
    # Return the final answer in the required format
    print(f"\n<<<{avg_layers:.4f}>>>")

solve_adsorption()