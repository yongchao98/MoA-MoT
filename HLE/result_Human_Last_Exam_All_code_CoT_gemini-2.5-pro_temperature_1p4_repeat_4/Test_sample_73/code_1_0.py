# Step 1: Define the Euler characteristic of the building blocks.

# The base space for our analysis is an annulus.
# An annulus is a disk with a smaller disk removed.
# chi(Disk) = 1.
# chi(Annulus) = chi(Disk) - chi(inner Disk removed) is not correct.
# A better way: an Annulus is topologically a cylinder.
# A cylinder has chi = 0.
chi_annulus = 0

# The configuration space is stated to be smooth, which means we must
# exclude self-intersecting configurations. These occur at 2 specific
# points in the annulus. So, we must "puncture" the annulus at these 2 points.
num_punctures = 2

# The Euler characteristic of a punctured surface is the original chi minus the number of punctures.
chi_base_space = chi_annulus - num_punctures

print(f"The base space X' is an annulus with {num_punctures} punctures.")
print(f"chi(Annulus) = {chi_annulus}")
print(f"chi(X') = chi(Annulus) - num_punctures = {chi_annulus} - {num_punctures} = {chi_base_space}")
print("-" * 20)

# Step 2: Construct the full configuration space M.
# M is a 2-sheeted cover of the base space X'.
# According to the theory of covering spaces, the Euler characteristic of the
# total space is the number of sheets times the Euler characteristic of the base space.
num_sheets = 2
chi_M = num_sheets * chi_base_space

print(f"The configuration space M is a {num_sheets}-sheeted cover of X'.")
print(f"chi(M) = num_sheets * chi(X') = {num_sheets} * {chi_base_space} = {chi_M}")
print("-" * 20)

# Step 3: Calculate the genus g from the Euler characteristic of M.
# The formula for a closed, orientable surface is chi(M) = 2 - 2*g.
# We need to solve for g.
# chi(M) = 2 - 2g
# 2g = 2 - chi(M)
# g = (2 - chi(M)) / 2

# Let's show the variables in the final equation
# 2 - 2g = chi_M
val_2 = 2
val_2g_coeff = 2

# 2g = 2 - chi_M
rhs = val_2 - chi_M
genus = rhs / val_2g_coeff

print("The genus g is found using the formula: 2 - 2*g = chi(M)")
print(f"Substituting the value for chi(M):")
print(f"{val_2} - {val_2g_coeff}*g = {chi_M}")
print(f"{val_2g_coeff}*g = {val_2} - ({chi_M})")
print(f"{val_2g_coeff}*g = {rhs}")
print(f"g = {rhs} / {val_2g_coeff}")
print(f"g = {int(genus)}")
