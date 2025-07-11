import math

def calculate_csfs():
    """
    Calculates the number of Configuration State Functions (CSFs) for a full CI
    calculation of CH2SiHF with a 6-31G** basis set.
    """

    # Step 1: Determine the total number of electrons (N)
    electrons = {
        'C': 6,
        'H': 1,
        'Si': 14,
        'F': 9
    }
    # Molecule is CH2SiHF, so 1*C, 3*H, 1*Si, 1*F
    num_atoms = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}
    N = sum(electrons[atom] * num for atom, num in num_atoms.items())

    print(f"Step 1: Calculate the total number of electrons (N)")
    print(f"Electrons per atom: C={electrons['C']}, H={electrons['H']}, Si={electrons['Si']}, F={electrons['F']}")
    print(f"Molecule CH2SiHF has 1 C, 3 H, 1 Si, 1 F atom(s).")
    print(f"Total electrons N = 1*{electrons['C']} + 3*{electrons['H']} + 1*{electrons['Si']} + 1*{electrons['F']} = {N}\n")

    # Step 2: Determine the total number of basis functions (K) for 6-31G**
    # 6-31G** basis functions per atom:
    # H: 2 (6-31G) + 3 (p-funcs from second *) = 5
    # C: 9 (6-31G) + 6 (d-funcs from first *) = 15
    # F: 9 (6-31G) + 6 (d-funcs from first *) = 15
    # Si: 13 (6-31G) + 6 (d-funcs from first *) = 19
    basis_functions = {
        'C': 15,
        'H': 5,
        'Si': 19,
        'F': 15
    }
    K = sum(basis_functions[atom] * num for atom, num in num_atoms.items())

    print(f"Step 2: Calculate the total number of basis functions (K)")
    print(f"Basis functions (orbitals) per atom for 6-31G**: C={basis_functions['C']}, H={basis_functions['H']}, Si={basis_functions['Si']}, F={basis_functions['F']}")
    print(f"Total basis functions K = 1*{basis_functions['C']} + 3*{basis_functions['H']} + 1*{basis_functions['Si']} + 1*{basis_functions['F']} = {K}\n")

    # Step 3: Apply the Weyl-Paldus formula for a singlet state (S=0)
    # Number of CSFs = (K+1) * [C(K+1, N/2)]^2
    # The ground state of this molecule is a singlet, so S=0.
    S = 0
    if N % 2 != 0:
        print("Error: The molecule has an odd number of electrons. The formula for S=0 is not applicable for the ground state.")
        return

    n_half = N // 2
    k_plus_1 = K + 1

    print("Step 3: Apply the Weyl-Paldus formula for a singlet (S=0) state")
    print("Formula: Number of CSFs = (K+1) * (C(K+1, N/2))^2")
    print(f"With N = {N} and K = {K}, we have N/2 = {n_half} and K+1 = {k_plus_1}.\n")

    # Step 4: Calculate the final value
    print("Step 4: Calculate the final result")
    print(f"The equation is: Number of CSFs = {k_plus_1} * (C({k_plus_1}, {n_half}))^2")
    
    # Calculate the binomial coefficient C(K+1, N/2)
    # Using math.comb for large number precision
    try:
        comb_val = math.comb(k_plus_1, n_half)
    except ValueError:
        print(f"Cannot calculate C({k_plus_1}, {n_half}) as n < k.")
        return
        
    # Calculate the final number of CSFs
    num_csfs = k_plus_1 * (comb_val ** 2)
    
    print(f"First, C({k_plus_1}, {n_half}) = {comb_val}")
    print(f"Then, the total number of CSFs = {k_plus_1} * ({comb_val})^2")
    print(f"Final Answer: The number of CSFs is {num_csfs}")

if __name__ == '__main__':
    calculate_csfs()