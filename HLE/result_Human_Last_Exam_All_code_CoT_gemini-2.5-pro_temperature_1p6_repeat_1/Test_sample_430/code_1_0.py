import math

def calculate_average_layers():
    """
    Calculates the average number of adsorbed layers per site using the BET model.
    
    Assumptions:
    1. The lateral interaction energy epsilon_l is negligible.
    2. The maximum number of layers k is large, so the infinite-layer BET model is applicable.
    3. The interaction energy for the 2nd layer and above (epsilon_v) is assumed to be 2 * epsilon_1
       to ensure a finite number of layers (mu < epsilon_v).
    """
    
    # Given parameters in units of k_B*T
    mu_kT = 0.15
    epsilon_1_kT = 0.1
    
    # Assumed parameter for the energy of subsequent layers (epsilon_v)
    epsilon_v_kT = 2 * epsilon_1_kT
    
    # Check if the condition for finite layers (mu < epsilon_v) is met
    if not (mu_kT < epsilon_v_kT):
        print("Warning: The condition mu < epsilon_v is not met. The number of layers would diverge.")
        return

    # Calculate the BET parameters x and c
    x = math.exp(mu_kT - epsilon_v_kT)
    c = math.exp(epsilon_v_kT - epsilon_1_kT)

    # Calculate the average number of layers using the BET equation
    # <h> = (c * x) / ((1 - x) * (1 - x + c * x))
    numerator = c * x
    denominator = (1 - x) * (1 - x + c * x)
    
    average_h = numerator / denominator

    # Output the final equation with all numbers plugged in
    print("Based on the BET model with the assumptions made, the calculation is:")
    print(f"<k> = (c * x) / ((1 - x) * (1 - x + c * x))")
    print(f"<k> = ({c:.4f} * {x:.4f}) / ((1 - {x:.4f}) * (1 - {x:.4f} + {c:.4f} * {x:.4f}))")
    print(f"<k> = ({numerator:.4f}) / (({1-x:.4f}) * ({1-x+c*x:.4f}))")
    print(f"<k> = {numerator:.4f} / {denominator:.4f}")
    print(f"The average number of adsorbed layers per site is: {average_h:.4f}")

calculate_average_layers()