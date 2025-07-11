import math

def solve():
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.
    """
    try:
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n < 1:
            print("Error: Number of qubits must be a positive integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer for n.")
        return

    # The number of destabilizer sets is given by the formula:
    # Total = (Number of phase choices) * (Number of Pauli string choices)
    # Number of phase choices = 4^n
    # Number of Pauli string choices = 2^(n*(n+1)/2)

    # Calculate the number of choices for the Pauli strings.
    # This is equivalent to counting the number of n x n symmetric binary matrices.
    # The exponent is n*(n+1)/2.
    exponent_pauli_strings = n * (n + 1) // 2
    
    # Using Python's arbitrary-precision integers to handle large numbers.
    num_pauli_strings = pow(2, exponent_pauli_strings)

    # Calculate the number of choices for the global phases.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    num_phase_choices = pow(4, n)

    # The total number is the product of the two.
    total_destabilizers = num_phase_choices * num_pauli_strings

    print(f"\nFor n = {n} qubits:")
    print("The final equation is: Total = (4^n) * (2^(n*(n+1)/2))")
    print("-" * 30)
    
    # Output each number in the final equation
    print(f"Number of phase choices (4^n): {num_phase_choices}")
    print(f"Number of Pauli string choices (2^({n}*({n}+1)/2)): {num_pauli_strings}")
    print(f"Total number of different destabilizer sets: {total_destabilizers}")

solve()
<<<4**n * 2**(n*(n+1)/2)>>>