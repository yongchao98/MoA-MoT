import math

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for theta.

    Args:
        d (int): The number of sensor nodes.
        F (float): The fidelity of the initial noisy GHZ state.
    """
    if not (isinstance(d, int) and d > 0):
        print("Error: d (number of nodes) must be a positive integer.")
        return
    if not (isinstance(F, (int, float)) and 0 <= F <= 1):
        print("Error: F (fidelity) must be a number between 0 and 1.")
        return

    # The Quantum Fisher Information (QFI) for the parameter theta is derived to be:
    # F_Q = 4 * d * (2*F - 1)^2
    qfi = 4 * d * (2 * F - 1)**2

    # The problem asks for the difference between 1 and the QFI.
    result = 1 - qfi

    # Print the final equation with all the numbers plugged in, as requested.
    # The numbers in the equation are 1, 4, d, 2, and F.
    print(f"For d = {d} and F = {F}:")
    print(f"The difference is calculated by the equation: 1 - 4 * d * (2 * F - 1)^2")
    print(f"Result = 1 - 4 * {d} * (2 * {F} - 1)**2 = {result}")

# --- Example Usage ---
# You can change the values of d and F below to explore different scenarios.

# d is the number of sensor nodes (must be a positive integer)
d = 10

# F is the fidelity of the initial state (must be between 0 and 1)
F = 0.95

calculate_qfi_difference(d, F)
