import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity of cyclopentanecarboxylic acid,
    the product of the Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    
    # The product is cyclopentanecarboxylic acid.
    # Its molecular graph consists of non-hydrogen atoms.
    # It has a 5-carbon ring, with a -COOH group attached to one carbon.

    # Step 1: Determine Na, the number of non-hydrogen atoms (vertices).
    # 5 carbons in the cyclopentane ring + 1 carbon in the carboxylic acid group
    # + 2 oxygens in the carboxylic acid group.
    N_a = 5 + 1 + 2

    # Step 2: Determine Nb, the number of bonds between non-hydrogen atoms (edges).
    # 5 C-C bonds in the ring.
    # 1 C-C bond connecting the ring to the carboxylic group.
    # 1 C=O and 1 C-O bond in the carboxylic group.
    N_b = 5 + 1 + 2

    # Step 3: Determine Nc, the number of cycles.
    # There is one cyclopentane ring.
    N_c = 1

    # Step 4: Determine the degree (δ) for each non-hydrogen atom and sum their squares.
    # Atom degrees (δ):
    # - Carbon on ring attached to COOH: bonded to 2 ring Cs, 1 carboxyl C -> δ=3
    # - 2 Carbons on ring adjacent to the first one: each bonded to 2 other ring Cs -> δ=2
    # - 2 other Carbons on ring: each bonded to 2 other ring Cs -> δ=2
    # - Carbon in COOH group: bonded to ring C, =O, and -OH -> δ=3
    # - Oxygen in =O group: bonded to C -> δ=1
    # - Oxygen in -OH group: bonded to C -> δ=1
    vertex_degrees = [3, 2, 2, 2, 2, 3, 1, 1]
    
    sum_sq_degrees = sum([d**2 for d in vertex_degrees])

    # Step 5: Calculate the Böttcher Molecular Complexity Index.
    # Formula: BMCI = (N_b^2 / N_a) * (N_c / N_a) * Σ(δ^2)
    complexity_index = (N_b**2 / N_a) * (N_c / N_a) * sum_sq_degrees

    print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.")
    print("Calculating its Böttcher Molecular Complexity:")
    print(f"Number of non-hydrogen atoms (Na): {N_a}")
    print(f"Number of bonds between non-hydrogen atoms (Nb): {N_b}")
    print(f"Number of cycles (Nc): {N_c}")
    print(f"Sum of squares of atomic degrees (Σ(δᵢᵛ)²): {sum_sq_degrees}")
    print("\nBöttcher Molecular Complexity Formula:")
    print("BMCI = (Nb² / Na) * (Nc / Na) * Σ(δᵢᵛ)²")
    print("\nFinal Calculation:")
    print(f"BMCI = ({N_b}² / {N_a}) * ({N_c} / {N_a}) * {sum_sq_degrees}")
    print(f"BMCI = ({N_b**2} / {N_a}) * ({N_c / N_a}) * {sum_sq_degrees}")
    print(f"BMCI = {N_b**2 / N_a} * {N_c / N_a} * {sum_sq_degrees}")
    print(f"BMCI = {complexity_index}")

calculate_bottcher_complexity()