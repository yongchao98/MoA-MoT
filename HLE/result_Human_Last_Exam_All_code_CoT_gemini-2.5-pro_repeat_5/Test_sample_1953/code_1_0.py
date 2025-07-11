import math

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI)
    for a distributed quantum sensing scenario.

    Args:
        d (int): The number of sensor nodes.
        F (float): The fidelity of the initial noisy GHZ state.
    """
    # Validate inputs
    if not isinstance(d, int) or d <= 0:
        print("Error: d (number of sensors) must be a positive integer.")
        return
    if not isinstance(F, (int, float)) or not (0 <= F <= 1):
        print("Error: F (fidelity) must be a number between 0 and 1.")
        return

    # The derived formula for the Quantum Fisher Information (QFI) is H_theta = 4 * d * (2*F - 1)^2.
    
    # --- Calculation Steps ---
    
    # Step 1: Calculate the term inside the parenthesis for clarity
    term_in_parenthesis = 2 * F - 1
    
    # Step 2: Calculate the Quantum Fisher Information (H_theta)
    qfi = 4 * d * (term_in_parenthesis**2)

    # Step 3: Calculate the final difference between 1 and the QFI
    result = 1 - qfi

    # --- Output the Results ---
    
    print("This script calculates 1 - QFI for the specified sensing scenario.")
    print(f"User-provided parameters: d = {d}, F = {F}")
    print("\n" + "="*40)
    print("Detailed Calculation:")
    print("="*40)

    # Output the formula for QFI
    print("The formula for the Quantum Fisher Information (QFI) is:")
    print("H_theta = 4 * d * (2*F - 1)^2")
    
    # Output the QFI calculation with substituted values
    print("\nSubstituting the values to find H_theta:")
    print(f"H_theta = 4 * {d} * (2 * {F} - 1)^2")
    print(f"H_theta = 4 * {d} * ({term_in_parenthesis})^2")
    print(f"H_theta = 4 * {d} * {term_in_parenthesis**2}")
    print(f"H_theta = {qfi}")

    # Output the final difference calculation
    print("\nCalculating the final result (1 - H_theta):")
    print(f"Result = 1 - {qfi}")
    print(f"Result = {result}")
    print("="*40)


# --- Example Execution ---
# You can change these values to explore different scenarios.

# Example 1: A 4-qubit system with high fidelity
d_example_1 = 4
F_example_1 = 0.9
calculate_qfi_difference(d_example_1, F_example_1)

print("\n\n")

# Example 2: A 10-qubit system with lower fidelity
d_example_2 = 10
F_example_2 = 0.75
calculate_qfi_difference(d_example_2, F_example_2)
