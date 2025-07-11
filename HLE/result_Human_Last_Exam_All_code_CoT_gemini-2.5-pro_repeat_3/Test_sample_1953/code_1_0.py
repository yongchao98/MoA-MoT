import numpy as np

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI)
    for a given number of sensor nodes d and fidelity F.

    Args:
        d (int): The number of sensor nodes (d).
        F (float): The fidelity of the initial state (F), between 0 and 1.
    """
    if not isinstance(d, int) or d <= 0:
        print("Error: d must be a positive integer.")
        return
    if not isinstance(F, (int, float)) or not (0 <= F <= 1):
        print("Error: F must be a number between 0 and 1.")
        return

    print(f"Given parameters:")
    print(f"Number of sensor nodes, d = {d}")
    print(f"Fidelity, F = {F}")
    print("-" * 20)
    
    # Calculate the Quantum Fisher Information (QFI)
    # QFI = 4 * d * (2*F - 1)^2
    term1 = 4 * d
    term2 = 2 * F - 1
    qfi = term1 * (term2 ** 2)

    # Calculate the final result: 1 - QFI
    result = 1 - qfi

    # Output the final equation with all the numbers
    print("The final result is calculated as: 1 - 4 * d * (2 * F - 1)^2")
    print(f"Result = 1 - 4 * {d} * (2 * {F} - 1)^2")
    print(f"Result = 1 - {term1} * ({term2})^2")
    print(f"Result = 1 - {term1} * {term2**2}")
    print(f"Result = 1 - {qfi}")
    print(f"Result = {result}")


if __name__ == '__main__':
    # Example values for d and F.
    # d is the number of qubits.
    # F is the fidelity.
    d_example = 4
    F_example = 0.9
    calculate_qfi_difference(d_example, F_example)
