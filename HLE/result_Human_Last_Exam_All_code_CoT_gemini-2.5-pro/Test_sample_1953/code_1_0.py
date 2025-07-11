import math

def calculate_fisher_information_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for a
    given number of sensor nodes (d) and initial state fidelity (F).

    Args:
        d (int): The number of sensor nodes (d >= 1).
        F (float): The fidelity of the initial noisy GHZ state (0 <= F <= 1).
    """
    if not isinstance(d, int) or d < 1:
        print("Error: 'd' must be an integer greater than or equal to 1.")
        return
    if not isinstance(F, (int, float)) or not (0 <= F <= 1):
        print("Error: 'F' must be a number between 0 and 1.")
        return

    print(f"Starting calculation for d = {d} and F = {F}")
    print("-" * 40)

    # The Quantum Fisher Information (I_theta) for the parameter theta is I_theta = 4 * d * (2*F - 1)^2
    # We need to output each number in the final equation.

    # Step 1: Calculate the term (2*F - 1)
    term1 = 2 * F - 1
    print(f"Step 1: Calculate the coherence term (2*F - 1).")
    print(f"   (2 * {F} - 1) = {term1:.4f}")

    # Step 2: Square the term
    term2 = term1 ** 2
    print(f"\nStep 2: Square the coherence term.")
    print(f"   ({term1:.4f})^2 = {term2:.4f}")

    # Step 3: Calculate the full QFI, I_theta
    I_theta = 4 * d * term2
    print(f"\nStep 3: Calculate the Quantum Fisher Information (I_theta).")
    print(f"   I_theta = 4 * d * (coherence term)^2")
    print(f"   I_theta = 4 * {d} * {term2:.4f} = {I_theta:.4f}")

    # Step 4: Calculate the final requested value: 1 - I_theta
    result = 1 - I_theta
    print(f"\nStep 4: Calculate the difference between 1 and the QFI.")
    print(f"   Result = 1 - I_theta")
    print(f"   Result = 1 - {I_theta:.4f} = {result:.4f}")

    print("-" * 40)
    print("Final equation with all numbers shown:")
    # Using f-string to print the equation with substituted numbers
    print(f"1 - (4 * {d} * (2 * {F} - 1)^2) = {result:.4f}")


# --- Example Usage ---
# You can change these values to explore different scenarios.
number_of_sensors = 5
fidelity = 0.95

calculate_fisher_information_difference(number_of_sensors, fidelity)