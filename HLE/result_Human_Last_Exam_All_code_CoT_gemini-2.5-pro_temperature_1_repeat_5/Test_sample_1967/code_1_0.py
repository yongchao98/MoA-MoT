def calculate_destabilizers():
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
        print("Error: Invalid input. Please enter a valid integer for n.")
        return

    # Calculate the number of choices for the Pauli part (number of n x n symmetric binary matrices)
    pauli_exponent = n * (n + 1) // 2
    
    # Calculate the number of choices for the phases (4 choices for each of the n generators)
    # In base 2, this is 2n
    phase_exponent_base2 = 2 * n
    
    # The total number of sets is the product, so we sum the exponents in base 2
    total_exponent = pauli_exponent + phase_exponent_base2
    
    # Calculate the final numbers
    num_pauli_choices = 2**pauli_exponent
    num_phase_choices = 4**n
    total_sets = 2**total_exponent

    print("-" * 50)
    print(f"For n = {n} qubits:")
    print("The total number of different destabilizer sets is the product of:")
    print(" 1. The number of choices for the Pauli operators.")
    print(" 2. The number of choices for the phases.")
    print()
    
    print("1. The number of Pauli operator choices equals the number of n x n symmetric binary matrices:")
    print(f"   2^(n*(n+1)/2) = 2^({n}*({n}+1)/2) = 2^{pauli_exponent} = {num_pauli_choices}")
    print()

    print("2. The number of phase choices for the n destabilizers is:")
    print(f"   4^n = 4^{n} = {num_phase_choices}")
    print()

    print("The total number of sets is the product:")
    print(f"   Total = (2^{pauli_exponent}) * (4^{n})")
    print(f"         = (2^{pauli_exponent}) * ((2^2)^{n})")
    print(f"         = 2^{pauli_exponent} * 2^{phase_exponent_base2}")
    print(f"         = 2^({pauli_exponent} + {phase_exponent_base2})")
    print(f"         = 2^{total_exponent}")
    print(f"         = {total_sets}")
    print("-" * 50)


if __name__ == '__main__':
    calculate_destabilizers()