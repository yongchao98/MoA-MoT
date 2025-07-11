import math

def calculate_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for theta.

    The QFI for theta is given by the formula: H(theta) = 4 * d * (2*F - 1)^2.

    Args:
        d (int): The number of sensor nodes (qubits).
        F (float): The fidelity of the noisy GHZ state with respect to the pure GHZ state.

    Returns:
        float: The value of 1 - H(theta).
    """
    # Check for valid inputs
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(F, (int, float)) or not (0 <= F <= 1):
        raise ValueError("F must be a float between 0 and 1.")

    # Define the numbers in the equation
    c_one = 1
    c_two = 2
    c_four = 4

    # Calculate the term (2*F - 1)
    coherence_term = c_two * F - c_one

    # Calculate the QFI
    qfi = c_four * d * (coherence_term ** c_two)

    # Calculate the final result
    result = c_one - qfi
    
    # Print the equation with the numbers plugged in
    print(f"Number of qubits d = {d}")
    print(f"Fidelity F = {F}")
    print("\nCalculating the Quantum Fisher Information, H(theta):")
    print(f"H(theta) = {c_four} * {d} * ({c_two} * {F} - {c_one})**{c_two}")
    print(f"H(theta) = {qfi}")
    
    print("\nCalculating the final result, 1 - H(theta):")
    print(f"Result = {c_one} - {qfi}")
    print(f"Result = {result}")
    
    return result

# --- Example Usage ---
# You can change these values to see the result for different scenarios.
d_value = 4  # Example: 4 qubits
F_value = 0.9 # Example: Fidelity of 90%

calculate_difference(d_value, F_value)
