import math

# Define parameters for the calculation
# Ne = 32 (electrons), Norb = 64 (orbitals for CH2SiHF in 6-31G**)
Ne = 32
Norb = 64

# For a singlet state (S=0), the number of alpha and beta electrons is n = Ne/2
n = Ne // 2

# The number of singlet CSFs is calculated using the formula:
# N_CSF(S=0) = C(Norb, n)^2 - C(Norb, n-1) * C(Norb, n+1)
# where C(n, k) is the binomial coefficient "n choose k".

# Calculate the necessary binomial coefficients
comb_n = math.comb(Norb, n)
comb_n_minus_1 = math.comb(Norb, n - 1)
comb_n_plus_1 = math.comb(Norb, n + 1)

# Calculate the two terms of the subtraction
term1 = comb_n * comb_n
term2 = comb_n_minus_1 * comb_n_plus_1

# Perform the final calculation
num_csfs = term1 - term2

# Print the final equation with all numbers and the result
print("The number of singlet (S=0) CSFs is given by the formula:")
print(f"Number of CSFs = C({Norb}, {n})^2 - C({Norb}, {n-1}) * C({Norb}, {n+1})")
print(f"Number of CSFs = {comb_n}^2 - {comb_n_minus_1} * {comb_n_plus_1}")
print(f"Number of CSFs = {term1} - {term2}")
print(f"Final number of CSFs = {num_csfs}")