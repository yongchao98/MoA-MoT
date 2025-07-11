import math

def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit
    stabilizer generator set {Z_1, ..., Z_n}.
    """
    # For demonstration, we'll use n=3 qubits. The logic is general for any n.
    n = 3

    # The stabilizer generator set is {S_1, ..., S_n} = {Z_1, ..., Z_n},
    # where Z_i is the Pauli Z operator on qubit i and Identity on others.
    # The destabilizer generator set {D_1, ..., D_n} must satisfy:
    # 1. [D_i, D_j] = 0         (destabilizers commute with each other)
    # 2. [S_i, D_j] = 0 for i!=j (commute with non-corresponding stabilizers)
    # 3. {S_i, D_i} = 0         (anti-commute with corresponding stabilizer)

    # Step 1: Determine the structure of each destabilizer D_j.
    # Let D_j be a tensor product of n Pauli operators.
    # From {S_j, D_j}=0, the Pauli operator on qubit j of D_j must anti-commute
    # with Z. Thus, it must be X or Y (2 choices).
    # From [S_i, D_j]=0 for i!=j, the Pauli operator on qubit i of D_j must
    # commute with Z. Thus, it must be I or Z (2 choices).

    # Step 2: Apply the mutual commutation condition [D_i, D_j]=0.
    # This condition imposes a symmetry constraint. For any pair (i, j),
    # the choice of Pauli (I or Z) on qubit j for D_i must be the same as
    # the choice of Pauli on qubit i for D_j.

    # Step 3: Count the number of valid combinations by parts.

    # Part A: Choices for the global phase of each generator.
    # Each of the n Pauli operators can be multiplied by a phase from {1, -1, i, -i}.
    # Number of phase choices = 4^n.
    phase_choices = 4**n

    # Part B: Choices for the "diagonal" Paulis.
    # For each generator D_i, the Pauli operator on qubit i must be X or Y.
    # Number of diagonal choices = 2^n.
    diagonal_choices = 2**n

    # Part C: Choices for the "off-diagonal" Paulis.
    # Due to the symmetry constraint, we only need to choose the Pauli (I or Z)
    # for each pair (i, j) where i < j. There are n*(n-1)/2 such pairs.
    num_off_diagonal_pairs = n * (n - 1) // 2
    off_diagonal_choices = 2**num_off_diagonal_pairs

    # The total number is the product of these independent choices.
    total_count = phase_choices * diagonal_choices * off_diagonal_choices

    # --- Output the results ---
    print(f"For an n={n} qubit system, the number of destabilizer sets is calculated as follows:")
    
    print("\n1. Number of choices for global phases (4 for each of n generators):")
    print(f"   4^{n} = {phase_choices}")
    
    print("\n2. Number of choices for diagonal Paulis (X or Y for each D_i on qubit i):")
    print(f"   2^{n} = {diagonal_choices}")

    print("\n3. Number of choices for symmetric off-diagonal Paulis (I or Z):")
    print(f"   The number of pairs (i,j) with i<j is {n}*({n}-1)/2 = {num_off_diagonal_pairs}")
    print(f"   2^{num_off_diagonal_pairs} = {off_diagonal_choices}")

    print("\n" + "="*50)
    print("The final equation for the total number of sets is the product:")
    print(f"Total = (Phase Choices) * (Diagonal Choices) * (Off-Diagonal Choices)")
    print(f"Total = 4^{n} * 2^{n} * 2^{num_off_diagonal_pairs}")
    print(f"Total = {phase_choices} * {diagonal_choices} * {off_diagonal_choices} = {total_count}")
    
    # The derived formula is 2^((n^2 + 5*n)/2).
    # Let's check for n=3: (3^2 + 5*3)/2 = (9+15)/2 = 12. 2^12 = 4096.

count_destabilizer_sets()
