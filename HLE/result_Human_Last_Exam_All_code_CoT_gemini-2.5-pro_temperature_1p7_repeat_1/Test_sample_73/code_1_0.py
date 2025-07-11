# Step 1: Define the Euler characteristics of the base topological spaces.
# The parameter space is a 2-torus, T^2.
chi_T2 = 0
# The branch curve B is a figure-eight, which has 2 vertices (nodes) and 4 edges.
V_B = 2
E_B = 4
chi_B = V_B - E_B
print(f"The Euler characteristic of the branch curve B is chi(B) = {V_B} - {E_B} = {chi_B}")

# The forbidden region R' is topologically a disk.
chi_R_prime = 1
print(f"The Euler characteristic of the forbidden region R' is chi(R') = {chi_R_prime}")

# Step 2: Calculate the Euler characteristic of the allowed region R.
# Using the formula: chi(T^2) = chi(R) + chi(R') - chi(B)
# 0 = chi(R) + 1 - (-2)
chi_R = chi_T2 - chi_R_prime + chi_B
print(f"The Euler characteristic of the allowed region R is chi(R) = {chi_T2} - {chi_R_prime} - ({chi_B}) = {chi_R}")

# Step 3: Calculate the Euler characteristic of the intermediate singular surface M_sing.
# M_sing is a double cover of R, branched over B.
# Formula: chi(M_sing) = 2 * chi(R) - chi(B)
chi_M_sing = 2 * chi_R - chi_B
print(f"The Euler characteristic of the singular surface M_sing is chi(M_sing) = 2 * {chi_R} - ({chi_B}) = {chi_M_sing}")

# Step 4: Calculate the Euler characteristic of the final smooth surface M.
# M is obtained by resolving the 2 singularities of M_sing.
# Resolving each singularity reduces the Euler characteristic by 1.
num_singularities = 2
chi_M = chi_M_sing - num_singularities
print(f"The Euler characteristic of the smooth surface M is chi(M) = {chi_M_sing} - {num_singularities} = {chi_M}")

# Step 5: Calculate the genus g from the Euler characteristic chi_M.
# Formula for an orientable surface: chi = 2 - 2g, so 2g = 2 - chi.
# We solve for the genus, g.
g = (2 - chi_M) / 2
print("\nThe final equation for the genus g is: 2 * g = 2 - chi(M)")
print(f"2 * g = 2 - ({chi_M})")
print(f"2 * g = {2 - chi_M}")
print(f"g = {(2 - chi_M)} / 2")
print(f"The genus of the configuration space is g = {int(g)}")
