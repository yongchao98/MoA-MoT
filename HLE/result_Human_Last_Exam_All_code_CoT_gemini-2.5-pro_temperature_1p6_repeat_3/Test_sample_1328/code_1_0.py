import math

def find_optimal_feedback_factor(total_power_2P, P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b' for a wireless channel.

    Args:
        total_power_2P (float): The total power budget (2P).
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the signal part of the second transmission.
        alpha (float): The weather-induced correlation between noise measurements.
    """
    print(f"Given System Parameters:")
    print(f"Total Power Budget (2P): {total_power_2P}")
    print(f"Power for first transmission (P1): {P1}")
    print(f"Power for second transmission signal (P2): {P2}")
    print(f"Noise correlation (alpha): {alpha}")
    print("-" * 25)

    # Validate inputs
    if P1 < 0 or P2 < 0:
        print("Error: Power values (P1, P2) cannot be negative.")
        return
    if not -1 < alpha < 1:
        print("Error: Correlation 'alpha' must be between -1 and 1 for a valid covariance matrix.")
        return

    # Step 1: Calculate remaining power for feedback adjustment term b^2
    print("Step 1: Calculate remaining power for the feedback term b^2.")
    print("The total power is tr(K_X) = P1 + P2 + b^2, which must be <= 2P.")
    print("To maximize information, we use the full power budget.")
    
    rem_power_b_squared = total_power_2P - P1 - P2
    
    print(f"Equation: b^2 = 2P - P1 - P2")
    print(f"b^2 = {total_power_2P} - {P1} - {P2} = {rem_power_b_squared}")
    print("-" * 25)

    if rem_power_b_squared < 0:
        print("Error: The sum of P1 and P2 exceeds the total power budget.")
        print("Cannot find a valid value for b.")
        return

    # Step 2: Calculate the maximum possible magnitude of b
    print("Step 2: Calculate the maximum magnitude of b, denoted as B.")
    max_b_magnitude = math.sqrt(rem_power_b_squared)
    print(f"Equation: B = sqrt(b^2)")
    print(f"B = sqrt({rem_power_b_squared}) = {max_b_magnitude}")
    print("-" * 25)

    # Step 3: Determine the optimal sign for b
    print("Step 3: Determine the optimal sign for b.")
    print("The objective is to maximize |K_{X+Z}|, which contains the term 2*P1*alpha*b.")
    print("To make this term positive and maximize the determinant, b should have the same sign as alpha.")

    if alpha >= 0:
        optimal_b = max_b_magnitude
        print(f"Since alpha ({alpha}) is non-negative, the optimal b should be positive.")
    else:
        optimal_b = -max_b_magnitude
        print(f"Since alpha ({alpha}) is negative, the optimal b should be negative.")
    
    print("-" * 25)
    print("Final Answer:")
    print(f"The optimal feedback adjustment factor b is: {optimal_b}")


# --- Example Usage ---
# You can change these values to see how the result changes.
if __name__ == '__main__':
    # Total power budget
    total_power_2P = 20.0
    # Power for the first transmission
    P1 = 8.0
    # Power for the independent part of the second transmission
    P2 = 10.0
    # Noise correlation
    alpha = -0.4  # Example with negative correlation
    
    find_optimal_feedback_factor(total_power_2P, P1, P2, alpha)

    print("\n" + "="*40 + "\n")

    # Another example with positive correlation
    alpha_pos = 0.7
    find_optimal_feedback_factor(total_power_2P, P1, P2, alpha_pos)
