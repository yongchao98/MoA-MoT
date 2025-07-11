import sys

def solve():
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.
    """
    try:
        # Prompt the user to enter the number of qubits, n.
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n <= 0:
            print("Error: Number of qubits (n) must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please enter an integer.", file=sys.stderr)
        return

    print(f"\nFor n = {n} qubits, we calculate the number of possible destabilizer sets.")
    print("The plan is to count:")
    print("1. The number of valid choices for the phase-less Pauli operators.")
    print("2. The number of ways to assign a phase to each operator.")
    
    # 1. Count the number of choices for the phase-less part of the operators.
    # As derived in the reasoning, this corresponds to the number of n x n
    # symmetric binary matrices. A symmetric matrix is defined by the elements
    # on and above its main diagonal.
    # Number of diagonal elements = n
    # Number of off-diagonal elements in the upper triangle = (n*n - n) / 2
    # Total number of elements to choose freely = n + (n*n - n)/2 = n*(n+1)/2
    num_free_params = n * (n + 1) // 2
    
    # Each of these parameters can be 0 or 1.
    num_pauli_sets = 2**num_free_params

    print("\nStep 1: Counting the phase-less Pauli operator sets.")
    print(f"The structure of the phase-less destabilizer generators is determined by an {n}x{n} symmetric binary matrix.")
    print(f"The number of free parameters in this matrix is n*(n+1)/2 = {n}*({n}+1)/2 = {num_free_params}.")
    print(f"Thus, the number of possible sets of phase-less Pauli operators is 2**{num_free_params} = {num_pauli_sets}.")

    # 2. Count the number of ways to assign phases.
    # Each of the n destabilizer generators is an element of the Pauli group,
    # so it can have a phase from {+1, -1, +i, -i}.
    num_phases = 4**n
    
    print("\nStep 2: Counting the phase choices.")
    print("Each of the n destabilizer generators can have a phase from {+1, -1, +i, -i}.")
    print(f"This gives 4 choices for each of the {n} generators.")
    print(f"Thus, the total number of ways to assign phases is 4**{n} = {num_phases}.")

    # 3. Calculate the total number of sets.
    # The total number is the product of the number of Pauli sets and the number of phase choices.
    total_sets = num_pauli_sets * num_phases

    print("\nFinal Calculation:")
    print("Total number of sets = (Number of phase-less sets) * (Number of phase choices)")
    print(f"Total number of sets = {num_pauli_sets} * {num_phases}")
    print(f"The final result is: {total_sets}")

if __name__ == '__main__':
    solve()