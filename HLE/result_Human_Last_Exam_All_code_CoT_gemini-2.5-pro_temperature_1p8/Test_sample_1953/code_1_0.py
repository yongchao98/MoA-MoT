import sys

def calculate_qfi_difference(d, F):
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI) for theta.

    The formula for the QFI is derived as I_Q(theta) = 4 * d * (2*F - 1)**2.
    The function computes 1 - I_Q(theta).

    Args:
        d (int): The number of sensor nodes.
        F (float): The fidelity of the initial state with the GHZ state |psi_+>.

    Returns:
        float: The value of 1 - I_Q(theta).
    """
    
    # Calculate the (2F-1)^2 term, which relates to the coherence of the mixed state.
    coherence_term = (2 * F - 1)**2
    
    # Calculate the Quantum Fisher Information for theta.
    I_Q = 4 * d * coherence_term
    
    # Calculate the final result: 1 - I_Q(theta).
    result = 1 - I_Q
    
    # As requested, output the numbers in the final equation.
    # The final equation is: 1 - (4 * d * (2*F - 1)^2) = result
    print("Derivation Steps:")
    print(f"1. Number of sensor nodes, d = {d}")
    print(f"2. Fidelity of the initial state, F = {F}")
    print(f"3. The QFI is given by I_Q(theta) = 4 * d * (2*F - 1)^2")
    print(f"   I_Q(theta) = 4 * {d} * (2*{F} - 1)^2 = 4 * {d} * {coherence_term} = {I_Q}")
    print("\nFinal Equation:")
    print(f"Difference = 1 - I_Q(theta)")
    print(f"Difference = 1 - {I_Q} = {result}")

    return result

if __name__ == "__main__":
    # Check for correct number of command-line arguments.
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <d> <F>")
        print("  <d>: The total number of sensor nodes (an integer).")
        print("  <F>: The fidelity of the initial state (a float between 0 and 1).")
        sys.exit(1)

    try:
        # Parse command-line arguments.
        d_val = int(sys.argv[1])
        F_val = float(sys.argv[2])
    except ValueError:
        print("Error: Invalid input.")
        print("  'd' must be an integer and 'F' must be a float.")
        sys.exit(1)

    # Validate the input values.
    if d_val <= 0:
        print("Error: The number of sensors 'd' must be a positive integer.")
        sys.exit(1)
    if not (0.0 <= F_val <= 1.0):
        print("Error: The fidelity 'F' must be between 0.0 and 1.0.")
        sys.exit(1)

    calculate_qfi_difference(d_val, F_val)
