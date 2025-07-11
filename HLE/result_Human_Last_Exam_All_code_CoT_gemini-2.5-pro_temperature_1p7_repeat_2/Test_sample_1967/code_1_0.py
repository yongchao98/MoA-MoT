import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.

    Returns:
        int: The total number of different destabilizer sets.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of qubits 'n' must be a positive integer.")
        return None

    # Part 1: Count the number of ways to choose the Pauli operators.
    # This corresponds to choosing an n x n symmetric binary matrix, which determines
    # whether the operator at a given position is (I or X) or (Z or Y).
    # The number of such matrices is 2^(n*(n+1)/2).
    exponent_pauli = n * (n + 1) // 2
    num_pauli_choices = 2**exponent_pauli

    # Part 2: Count the number of ways to choose the global phases.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    # This gives 4^n possible combinations of phases.
    num_phase_choices = 4**n

    # The total number is the product of the two parts.
    total_count = num_pauli_choices * num_phase_choices

    return total_count, exponent_pauli, num_pauli_choices, num_phase_choices

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # Set the number of qubits here
    n = 3

    result = count_destabilizer_sets(n)

    if result:
        total_count, exponent_pauli, num_pauli_choices, num_phase_choices = result
        print(f"For an n = {n} qubit system with stabilizers {{Z_1, ..., Z_n}}:")
        print("\nStep 1: Calculate the number of choices for the Pauli operators.")
        print(f"This is given by the formula 2^(n*(n+1)/2).")
        print(f"Number of Pauli choices = 2^({n}*({n}+1)/2) = 2^{exponent_pauli} = {num_pauli_choices}")

        print("\nStep 2: Calculate the number of choices for the global phases.")
        print(f"This is given by the formula 4^n.")
        print(f"Number of phase choices = 4^{n} = {num_phase_choices}")

        print("\nStep 3: Calculate the total number of different destabilizer sets.")
        print("Total number = (Number of Pauli choices) * (Number of phase choices)")
        print(f"Total number = {num_pauli_choices} * {num_phase_choices} = {total_count}")

if __name__ == "__main__":
    main()