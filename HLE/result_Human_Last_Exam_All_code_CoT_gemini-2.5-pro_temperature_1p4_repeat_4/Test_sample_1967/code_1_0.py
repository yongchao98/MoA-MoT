def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    try:
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n <= 0:
            print("Error: The number of qubits (n) must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Part 1: Calculate the number of phase choices.
    # Each of the n destabilizer generators can have a phase from {1, -1, i, -i}.
    num_phases = 4**n

    # Part 2: Calculate the number of choices for the Pauli operators.
    # This is equal to the number of n x n symmetric binary matrices.
    # The number of freely chosen elements in such a matrix is n(n+1)/2.
    exponent = n * (n + 1) // 2
    num_pauli_choices = 2**exponent

    # Part 3: The total number of sets is the product of the two parts.
    total_sets = num_phases * num_pauli_choices

    # Output the explanation and the final equation.
    print(f"\nFor n = {n}, the calculation is as follows:")
    print(f"1. Number of phase choices for the n generators: 4^{n} = {num_phases}")
    print(f"2. Number of choices for the Pauli operator structures: 2^(n*(n+1)/2) = 2^{exponent} = {num_pauli_choices}")
    print("\nThe total number of different destabilizer sets is the product of these two numbers:")
    print(f"{num_phases} * {num_pauli_choices} = {total_sets}")

if __name__ == '__main__':
    count_destabilizer_sets()