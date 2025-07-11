import math

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for theta.

    The problem describes a distributed quantum sensing scenario with d sensors.
    The parameter to estimate is theta = sum(x_i) / sqrt(d).
    The initial state is a noisy GHZ state with fidelity F with respect to the
    (|0...0> + |1...1>)/sqrt(2) state.

    The derived formula for the QFI for theta is: H(theta) = 4 * d * (2*F^2 - 1)^2.
    This function calculates 1 - H(theta).
    
    Args:
        d (int): The number of sensor nodes. Must be a positive integer.
        F (float): The fidelity of the noisy GHZ state. Must be between 0.0 and 1.0.

    Returns:
        float: The value of 1 - QFI.
    """
    if not isinstance(d, int) or d <= 0:
        raise ValueError("d must be a positive integer.")
    if not isinstance(F, (int, float)) or not (0.0 <= F <= 1.0):
        raise ValueError("Fidelity F must be a number between 0 and 1.")

    # The QFI H(theta) is given by the formula 4 * d * (2*F^2 - 1)^2
    term_1 = 4
    term_2 = d
    # The term (2F^2 - 1) is related to the purity of the state in the relevant basis.
    fidelity_term_sq = (2 * F**2 - 1)**2
    
    qfi = term_1 * term_2 * fidelity_term_sq
    
    result = 1 - qfi
    
    print("This script calculates the difference between 1 and the Quantum Fisher Information (QFI).")
    print(f"For the given parameters: d = {d}, F = {F}")
    print("-" * 30)
    print("The final equation for the result is: 1 - (4 * d * (2 * F^2 - 1)^2)")
    print("\nBreaking down the calculation:")
    print(f"Quantum Fisher Information (QFI) = {term_1} * {term_2} * (2 * {F}^2 - 1)^2")
    print(f"QFI = {term_1} * {term_2} * {fidelity_term_sq:.4f}")
    print(f"QFI = {qfi:.4f}")
    print("\nFinal Result:")
    # The instruction "output each number in the final equation" is followed here
    print(f"Difference = 1 - QFI = {1} - {qfi:.4f} = {result:.4f}")
    
    return result

if __name__ == '__main__':
    # Example values for d and F.
    # d is the number of sensor nodes.
    # F is the fidelity of the initial state.
    num_sensors = 10
    fidelity = 0.95
    
    calculate_qfi_difference(num_sensors, fidelity)
