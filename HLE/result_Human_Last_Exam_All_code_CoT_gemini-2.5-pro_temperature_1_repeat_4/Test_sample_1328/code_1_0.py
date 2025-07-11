import numpy as np

def calculate_optimal_feedback(P1, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    The mutual information is maximized when the determinant of the covariance matrix
    of the received signal, |K_{X+Z}|, is maximized. The derivation shows that
    |K_{X+Z}| is a quadratic function of 'b':
    |K_{X+Z}| = -b^2 - (2 * P1 * alpha) * b + constant_terms

    To maximize this function, we find the vertex of the parabola by taking the
    derivative with respect to 'b' and setting it to zero:
    d/db(|K_{X+Z}|) = -2*b - 2*P1*alpha = 0
    This gives the optimal b = -P1 * alpha.

    This solution is valid under the condition that the power of the information
    signal in the second transmission (P2') is non-negative, which requires
    P2 >= (P1^2 * alpha^2).

    Args:
        P1 (float): Power used in the first transmission.
        alpha (float): Weather-induced correlation between noise measurements.

    Returns:
        float: The optimal feedback adjustment factor 'b'.
    """
    if not (P1 >= 0):
        raise ValueError("Power P1 must be non-negative.")
    if not (-1 < alpha < 1):
        raise ValueError("Correlation alpha must be between -1 and 1.")

    b_optimal = -P1 * alpha
    return b_optimal

# --- Example Usage ---
# Given parameters for the system
# Power in the first transmission
P1 = 10.0
# Correlation factor
alpha = 0.5

# Calculate the optimal feedback factor
b = calculate_optimal_feedback(P1, alpha)

# Print the final result showing the equation
print("Objective: Find the optimal feedback adjustment factor 'b'.")
print("The formula for the optimal 'b' is derived as: b = -P1 * alpha")
print("\nGiven the values:")
print(f"P1 = {P1}")
print(f"alpha = {alpha}")
print("\nThe calculation is:")
# The final print statement shows the numbers in the equation
print(f"b = -({P1}) * ({alpha}) = {b}")
