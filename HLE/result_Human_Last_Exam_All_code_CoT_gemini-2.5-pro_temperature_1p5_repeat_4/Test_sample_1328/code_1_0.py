import math

def calculate_optimal_b(alpha, P, P1):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        alpha (float): The weather-induced correlation between noise measurements.
        P (float): Half of the total power budget (Total budget = 2P).
        P1 (float): Power used in the first transmission.
    """
    
    print(f"Given parameters: alpha = {alpha}, P = {P}, P1 = {P1}\n")

    # Validate inputs
    if not -1 < alpha < 1:
        print("Error: alpha must be between -1 and 1 for a valid covariance matrix.")
        return
    
    power_for_b_squared = 2 * P - P1
    if power_for_b_squared < 0:
        print("Error: Power in first transmission (P1) exceeds the total power budget (2P).")
        return

    # Determine the sign of b. It should be opposite to the sign of alpha
    # to maximize the distance |b - alpha|.
    # If alpha is 0, the choice of sign is arbitrary; we'll choose positive.
    if alpha > 0:
        sign_b = -1.0
    else: # This covers alpha <= 0
        sign_b = 1.0

    # Calculate the magnitude of b
    magnitude_b = math.sqrt(power_for_b_squared)

    # Calculate the optimal b
    optimal_b = sign_b * magnitude_b

    # --- Output the calculation steps ---
    print("Objective: Maximize Mutual Information, which simplifies to maximizing (b - alpha)^2.")
    print(f"Constraint: Power constraint P1 + b^2 <= 2*P leads to b^2 <= 2*{P} - {P1}, so b^2 <= {power_for_b_squared}.")
    print(f"This means |b| <= sqrt({power_for_b_squared}) = {magnitude_b:.4f}.")
    print("\nTo maximize (b - alpha)^2, we choose b at the boundary of its allowed range and with a sign opposite to alpha.")
    
    print("\nFinal Calculation:")
    print(f"b = -sgn(alpha) * sqrt(2*P - P1)")
    print(f"b = {-1 if sign_b > 0 else 1 if sign_b < 0 else 'N/A'} * sqrt(2*{P} - {P1})") # sgn() definition is opposite of sign_b
    print(f"b = {sign_b} * {magnitude_b:.4f}")
    print(f"Optimal b = {optimal_b:.4f}")

# Example Usage with sample values
# Let's assume a correlation of 0.5, a total power budget of 20 (so P=10),
# and a first transmission power of 4.
alpha_val = 0.5
P_val = 10
P1_val = 4

calculate_optimal_b(alpha_val, P_val, P1_val)
