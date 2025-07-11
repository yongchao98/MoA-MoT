import math

def find_optimal_b(P, P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P (float): Half of the total power budget (budget = 2P).
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission message signal.
        alpha (float): Weather-induced correlation between noise measurements.
    """
    # The final equation for the optimal b is derived from maximizing the mutual information
    # under the power constraint tr(K_X) <= 2P.
    final_equation = "b = sign(alpha) * sqrt(2*P - P1 - P2)"

    # Check if the power allocation is physically possible
    if P1 < 0 or P2 < 0:
        print("Error: Power values P1 and P2 cannot be negative.")
        return
        
    power_budget = 2 * P
    power_consumed_by_signals = P1 + P2
    
    if power_consumed_by_signals > power_budget:
        print(f"Error: Power allocation P1 + P2 = {power_consumed_by_signals} exceeds the total budget 2P = {power_budget}.")
        return

    # Calculate the value under the square root from the power constraint
    # b^2 <= 2P - P1 - P2
    val_under_sqrt = power_budget - power_consumed_by_signals
    
    # The maximum magnitude of b is B
    B = math.sqrt(val_under_sqrt)

    # Determine the sign of b based on the sign of alpha
    if alpha >= 0:
        b_optimal = B
        sign_alpha_str = "1 (since alpha >= 0)"
    else:
        b_optimal = -B
        sign_alpha_str = "-1 (since alpha < 0)"

    print("The optimal feedback adjustment factor 'b' is calculated to maximize channel capacity.")
    print("The final equation for optimal b is:")
    print(f"  {final_equation}\n")
    print("Substituting the given parameter values:")
    print(f"  b = sign({alpha}) * sqrt(2*{P} - {P1} - {P2})")
    print(f"  b = {sign_alpha_str} * sqrt({val_under_sqrt})")
    print(f"  b = {sign_alpha_str} * {B:.4f}")
    print("\n" + "="*40)
    print(f"The final result for the optimal b is: {b_optimal:.4f}")
    print("="*40)

# --- Example Usage ---
# Given system parameters
P = 10.0   # Total power budget = 2P = 20
P1 = 8.0   # Power in first transmission
P2 = 5.0   # Power in second transmission's message
alpha = 0.6 # Positive noise correlation

find_optimal_b(P, P1, P2, alpha)

# --- Second Example ---
print("\n--- Another Example ---")
alpha_negative = -0.4 # Negative noise correlation
find_optimal_b(P=25.0, P1=15.0, P2=20.0, alpha=alpha_negative)