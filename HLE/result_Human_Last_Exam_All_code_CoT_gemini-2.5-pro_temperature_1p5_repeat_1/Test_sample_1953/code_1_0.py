import math

def solve_qfi_problem(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for a
    given distributed quantum sensing scenario.

    Args:
        d (int): The total number of sensor nodes.
        F (float): The fidelity of the noisy d-qubit GHZ state.
                   It must be between 0 and 1.
    """
    if not (isinstance(d, int) and d > 0):
        print("Error: 'd' must be a positive integer.")
        return
    if not (isinstance(F, (int, float)) and 0 <= F <= 1):
        print("Error: 'F' must be a number between 0 and 1.")
        return

    print(f"Calculating for d = {d} and F = {F}\n")

    # The Quantum Fisher Information (QFI) for this problem is given by the formula:
    # H(theta) = 4 * d * (2 * F^2 - 1)^2
    
    # Let's break down the calculation.
    # Term 1: 4
    term1 = 4
    
    # Term 2: d
    term2 = d
    
    # Term 3: F-squared
    F_squared = F**2
    
    # Term 4: (2 * F^2 - 1)
    fidelity_factor_base = 2 * F_squared - 1
    
    # Term 5: (2 * F^2 - 1)^2
    fidelity_factor_squared = fidelity_factor_base**2
    
    # Calculate the QFI
    qfi = term1 * term2 * fidelity_factor_squared
    
    # Calculate the final result
    result = 1 - qfi

    # Output the detailed equation
    print("The Quantum Fisher Information H(theta) is calculated as: 4 * d * (2 * F^2 - 1)^2")
    print(f"H(theta) = {term1} * {d} * (2 * {F:.2f}^2 - 1)^2")
    print(f"H(theta) = {term1} * {d} * (2 * {F_squared:.4f} - 1)^2")
    print(f"H(theta) = {term1} * {d} * ({fidelity_factor_base:.4f})^2")
    print(f"H(theta) = {term1} * {d} * {fidelity_factor_squared:.4f}")
    print(f"H(theta) = {qfi:.4f}\n")
    
    print("The difference between 1 and the QFI is:")
    # Print the final equation with all numbers
    print(f"1 - H(theta) = 1 - {qfi:.4f}")
    print(f"Final Result = {result:.4f}")


# Main execution block
if __name__ == "__main__":
    # Assigning values for d and F as they were not specified in the problem.
    # d is the number of sensor nodes.
    # F is the fidelity.
    d_value = 4
    F_value = 0.9
    
    solve_qfi_problem(d_value, F_value)