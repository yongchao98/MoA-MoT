import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Step 1: Count the number of possible Pauli structures.
    # This is determined by the number of n x n symmetric binary matrices.
    # The number of independent elements in such a matrix is n*(n+1)/2.
    pauli_structure_exponent = n * (n + 1) // 2
    
    # Using python's arbitrary-precision integers for large results
    num_pauli_structures = 2**pauli_structure_exponent

    # Step 2: Count the number of phase combinations.
    # Each of the n destabilizers can have a phase from {1, -1, i, -i}.
    phase_choices_per_destabilizer = 4
    num_phase_combinations = phase_choices_per_destabilizer**n

    # Step 3: Calculate the total number.
    total_sets = num_pauli_structures * num_phase_combinations

    # Step 4: Express the entire calculation and result clearly.
    # The final formula is 2^(n(n+1)/2) * 4^n = 2^((n^2+5n)/2)
    final_exponent = (n**2 + 5 * n) // 2

    print(f"For an n={n} qubit system with stabilizer generators {{Z_1, ..., Z_n}}:\n")
    
    print("1. Counting the possible Pauli operator structures (ignoring phases):")
    print(f"   The number of valid structures corresponds to the number of {n}x{n} symmetric binary matrices.")
    print(f"   This is given by the formula 2^(n*(n+1)/2).")
    print(f"   Calculation: 2^({n}*({n}+1)/2) = 2^{pauli_structure_exponent}")
    if n <= 5: # Only print smaller results
         print(f"   Result: {num_pauli_structures}\n")
    else:
         print("\n")


    print("2. Counting the phase combinations for each structure:")
    print(f"   Each of the {n} destabilizer operators can have one of 4 phases (±1, ±i).")
    print(f"   This is given by the formula 4^n.")
    print(f"   Calculation: {phase_choices_per_destabilizer}^{n}")
    if n <= 8: # Only print smaller results
        print(f"   Result: {num_phase_combinations}\n")
    else:
        print("\n")
        
    print("3. Calculating the total number of different destabilizer sets:")
    print("   Total = (Number of Pauli Structures) * (Number of Phase Combinations)")
    print(f"   Total = (2^{pauli_structure_exponent}) * ({phase_choices_per_destabilizer}^{n})")
    
    print("\nThe final expression can be simplified to a single power of 2:")
    print("   Total = 2^((n^2 + 5*n)/2)")
    print(f"   Calculation: 2^(({n}^2 + 5*{n})/2) = 2^{final_exponent}")
    
    print("-" * 30)
    print(f"The final answer for n={n} is: {total_sets}")
    print("-" * 30)


# --- User Execution ---
# You can change the value of n to see the result for a different number of qubits.
# For example, let's calculate it for n=3 qubits.
n = 3
count_destabilizer_sets(n)