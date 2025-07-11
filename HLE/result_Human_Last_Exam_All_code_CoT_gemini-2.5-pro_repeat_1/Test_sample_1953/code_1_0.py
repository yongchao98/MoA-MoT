import sys

def solve_quantum_sensing():
    """
    Calculates the difference between 1 and the Quantum Fisher Information (QFI)
    for a parameter theta in a distributed quantum sensing scenario.

    The QFI for theta is derived to be H(theta) = 4 * d * (2*F - 1)^2,
    where d is the number of sensor nodes and F is the fidelity of the
    initial noisy GHZ state.

    This script prompts the user for d and F, calculates 1 - H(theta),
    and prints the result along with the full equation.
    """
    try:
        # Prompt user for the number of sensor nodes, d
        d_str = input("Enter the number of sensor nodes (d): ")
        d = int(d_str)
        if d <= 0:
            print("Error: The number of nodes 'd' must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input for 'd'. Please enter an integer.", file=sys.stderr)
        return

    try:
        # Prompt user for the fidelity, F
        f_str = input("Enter the fidelity of the initial state (F): ")
        F = float(f_str)
        if not (0 <= F <= 1):
            print("Error: The fidelity 'F' must be a value between 0 and 1 (inclusive).", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input for 'F'. Please enter a number.", file=sys.stderr)
        return

    # The QFI H(theta) is 4 * d * (2*F - 1)^2
    # First, calculate the term (2*F - 1)^2
    fidelity_term = (2 * F - 1)**2
    
    # Second, calculate the Quantum Fisher Information H(theta)
    qfi = 4 * d * fidelity_term
    
    # Finally, calculate the required difference: 1 - H(theta)
    result = 1 - qfi
    
    # Output the final equation with all the numbers, as requested.
    print("\nBased on the derivation, the Quantum Fisher Information H(theta) is given by the formula: 4 * d * (2*F - 1)^2")
    print(f"For d = {d} and F = {F}, the QFI is: 4 * {d} * (2 * {F} - 1)^2 = {qfi}")
    print("\nThe difference between 1 and the QFI is:")
    print(f"1 - H(theta) = 1 - {qfi} = {result}")

if __name__ == "__main__":
    solve_quantum_sensing()
