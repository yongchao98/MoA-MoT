import math

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for theta.

    Args:
        d (int): The total number of sensor nodes.
        F (float): The fidelity of the initial state with respect to the GHZ+ state.
                   Should be between 0 and 1.
    """
    if not isinstance(d, int) or d <= 0:
        print("Error: d (number of sensors) must be a positive integer.")
        return
    if not isinstance(F, (int, float)) or not (0 <= F <= 1):
        print("Error: F (fidelity) must be a number between 0 and 1.")
        return

    # The formula for the Quantum Fisher Information for theta is H(theta) = 4*d*(2*F-1)^2
    qfi = 4 * d * (2 * F - 1)**2
    
    # The final result is the difference between 1 and the QFI.
    difference = 1 - qfi

    print(f"Given parameters:")
    print(f"Number of sensor nodes (d) = {d}")
    print(f"Fidelity (F) = {F}")
    print("-" * 30)
    
    print("Calculation Steps:")
    print(f"1. Quantum Fisher Information (QFI) for theta is calculated using the formula:")
    print(f"   QFI = 4 * d * (2 * F - 1)^2")
    print(f"   QFI = 4 * {d} * (2 * {F} - 1)^2 = {qfi:.4f}")
    
    print("\n2. The difference is calculated as 1 - QFI:")
    print(f"   Difference = 1 - {qfi:.4f} = {difference:.4f}")
    
    print("-" * 30)
    print("Final equation with values substituted:")
    print(f"Result = 1 - (4 * {d} * (2 * {F} - 1)^2) = {difference:.4f}")

if __name__ == '__main__':
    # Example values for d and F.
    # You can change these values to see the results for different scenarios.
    num_sensors = 10
    fidelity = 0.95
    
    calculate_qfi_difference(num_sensors, fidelity)
    
    # Example with F=1 (perfect GHZ+ state)
    # print("\n" + "="*30 + "\n")
    # calculate_qfi_difference(d=10, F=1.0)
    
    # Example with F=0.5 (maximally mixed state in the GHZ subspace)
    # print("\n" + "="*30 + "\n")
    # calculate_qfi_difference(d=10, F=0.5)