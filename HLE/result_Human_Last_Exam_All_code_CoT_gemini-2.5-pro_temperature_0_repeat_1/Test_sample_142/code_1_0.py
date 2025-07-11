import math

# Step 1 & 2: Determine the number of electrons (N) and spin (S)
# Molecule: CH2SiHF
# Electrons from C: 6
# Electrons from H: 1 (x3)
# Electrons from Si: 14
# Electrons from F: 9
num_electrons = 6 + (2 * 1) + 14 + 1 + 9
# For a singlet ground state, the total spin quantum number S is 0.
spin_S = 0

# Step 3: Determine the number of basis functions (K) for the 6-31G** basis set
# C in 6-31G**: 14
# H in 6-31G**: 5
# Si in 6-31G**: 18
# F in 6-31G**: 14
# Total K = 1*C + 3*H + 1*Si + 1*F
num_orbitals_K = (1 * 14) + (3 * 5) + (1 * 18) + (1 * 14)

# Step 4: Calculate the number of CSFs using the genealogical formula for S=0
# No. CSFs = C(K, N/2)^2 - C(K, N/2 - 1) * C(K, N/2 + 1)
N_half = num_electrons // 2

# Calculate the individual combination terms
# C(K, N/2)
comb_k = math.comb(num_orbitals_K, N_half)
# C(K, N/2 - 1)
comb_k_minus_1 = math.comb(num_orbitals_K, N_half - 1)
# C(K, N/2 + 1)
comb_k_plus_1 = math.comb(num_orbitals_K, N_half + 1)

# Calculate the final number of CSFs
num_csfs = comb_k**2 - comb_k_minus_1 * comb_k_plus_1

# Print the results and the final equation
print("Calculation for the number of CSFs in a Full CI for CH2SiHF with 6-31G** basis set")
print("-" * 80)
print(f"Number of electrons (N): {num_electrons}")
print(f"Number of spatial orbitals (K): {num_orbitals_K}")
print(f"Total spin (S): {spin_S} (Singlet State)")
print(f"Value of N/2: {N_half}")
print("-" * 80)
print("Using the formula: No. CSFs = C(K, N/2)^2 - C(K, N/2 - 1) * C(K, N/2 + 1)")
print(f"Plugging in the values: No. CSFs = C({num_orbitals_K}, {N_half})^2 - C({num_orbitals_K}, {N_half - 1}) * C({num_orbitals_K}, {N_half + 1})")
print("-" * 80)
print(f"Intermediate values:")
print(f"C({num_orbitals_K}, {N_half}) = {comb_k}")
print(f"C({num_orbitals_K}, {N_half - 1}) = {comb_k_minus_1}")
print(f"C({num_orbitals_K}, {N_half + 1}) = {comb_k_plus_1}")
print("-" * 80)
print("Final Equation:")
print(f"Number of CSFs = {comb_k}^2 - {comb_k_minus_1} * {comb_k_plus_1}")
print(f"Number of CSFs = {num_csfs}")
