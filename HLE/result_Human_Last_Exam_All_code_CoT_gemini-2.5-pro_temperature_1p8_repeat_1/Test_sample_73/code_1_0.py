# Step 1: Define the properties of the parameter space S.
# The parameter space of the two main arms is a 2-torus.
genus_of_parameter_space = 1

# The boundary of the valid region S on the torus consists of 3 disjoint curves.
# This means S is a torus with 3 holes.
num_boundary_curves = 3

print(f"The parameter space S is a surface of genus {genus_of_parameter_space} with {num_boundary_curves} boundary components.")

# Step 2: Calculate the Euler characteristic of S.
# The formula for a surface with genus g_s and h boundaries is chi = 2 - 2*g_s - h.
chi_S = 2 - 2 * genus_of_parameter_space - num_boundary_curves
print(f"The Euler characteristic of S is: chi(S) = 2 - 2*{genus_of_parameter_space} - {num_boundary_curves} = {chi_S}")

# Step 3: Calculate the Euler characteristic of the configuration space M.
# The configuration space M is the double of S. Its Euler characteristic is 2 * chi(S).
chi_M = 2 * chi_S
print(f"The configuration space M is the double of S, so its Euler characteristic is: chi(M) = 2 * {chi_S} = {chi_M}")

# Step 4: Calculate the genus of M.
# For a closed orientable surface, chi = 2 - 2g. So, g = (2 - chi) / 2.
genus_M = (2 - chi_M) / 2
print(f"The genus of the configuration space M is calculated using g = (2 - chi(M)) / 2.")
print(f"g = (2 - ({chi_M})) / 2 = {int(genus_M)}")
