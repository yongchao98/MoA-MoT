import math

def calculate_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit 
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    if n <= 0:
        print("Number of qubits (n) must be a positive integer.")
        return

    # Part 1: Number of choices for the Z-part structure.
    # This corresponds to the number of n x n symmetric binary matrices
    # with an undefined diagonal, which is 2^(n*(n-1)/2).
    num_z_structures = 2**(n * (n - 1) // 2)

    # Part 2: Number of choices for the X/Y part.
    # For each of the n destabilizers, the main operator can be X or Y.
    num_xy_choices = 2**n

    # Part 3: Number of phase choices.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    num_phase_choices = 4**n
    
    # The final equation is the product of these three parts.
    total_sets = num_z_structures * num_xy_choices * num_phase_choices

    print(f"For an {n}-qubit system with stabilizers {{Z_1, ..., Z_{n}}}:")
    print("\nThe calculation is broken down into three parts:")
    
    # Outputting each number in the final equation
    print("\n1. Number of valid Z-part structures:")
    print(f"   Equation: 2^(n * (n - 1) / 2)")
    print(f"   Calculation: 2^({n} * {n-1} / 2) = 2^{n*(n-1)//2} = {num_z_structures}")

    print("\n2. Number of X/Y-part choices:")
    print(f"   Equation: 2^n")
    print(f"   Calculation: 2^{n} = {num_xy_choices}")
    
    print("\n3. Number of phase choices:")
    print(f"   Equation: 4^n")
    print(f"   Calculation: 4^{n} = {num_phase_choices}")
    
    print("\nTotal number of different destabilizer sets is the product of these parts:")
    print(f"   Final Equation: (Z-structures) * (X/Y choices) * (Phase choices)")
    print(f"   Calculation: {num_z_structures} * {num_xy_choices} * {num_phase_choices} = {total_sets}")
    
    # As a check, we can use the combined formula: 2^((n^2 + 5n)/2)
    exponent = (n**2 + 5 * n) // 2
    check_total = 2**exponent
    print(f"\nSimplified formula check: 2^((n^2 + 5n) / 2) = 2^(({n**2} + {5*n}) / 2) = 2^{exponent} = {check_total}")

# --- User Input ---
# You can change the value of n here to see the result for a different number of qubits.
n_qubits = 3
calculate_destabilizer_sets(n_qubits)

# Final answer is the calculated number
final_answer = 2**((n_qubits**2 + 5 * n_qubits) // 2)
# To match the output format, converting the integer to a string.
print(f'<<<{final_answer}>>>')