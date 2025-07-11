# Final calculation based on the derivation.

# The Euler characteristic of the region R (a punctured torus).
chi_R = -1

# The Euler characteristic of the boundary of R (a circle).
chi_boundary_R = 0

# The Euler characteristic of the configuration space M is given by
# the Riemann-Hurwitz formula for this type of double cover.
chi_M = 2 * chi_R - chi_boundary_R

# The genus 'g' of a closed orientable surface is related to its
# Euler characteristic by chi = 2 - 2g.
# So, g = (2 - chi_M) / 2.
g = (2 - chi_M) / 2

# We print the equation that leads to the final answer
# The equation is: chi_M = 2 - 2g
# Substituting the calculated value of chi_M:
# -2 = 2 - 2g
# This simplifies to 2g = 2 - (-2) = 4
# And finally g = 4 / 2 = 2

print("The Euler characteristic of the configuration space M is Chi(M) = 2 * Chi(R) - Chi(∂R).")
print(f"Chi(R) = {chi_R}")
print(f"Chi(∂R) = {chi_boundary_R}")
print(f"So, Chi(M) = 2 * {chi_R} - {chi_boundary_R} = {chi_M}")
print("\nThe genus 'g' is found using the formula Chi(M) = 2 - 2g.")
print(f"{chi_M} = 2 - 2g")
print(f"2g = 2 - ({chi_M})")
print(f"2g = {2 - chi_M}")
print(f"g = {2 - chi_M} / 2")
print(f"g = {int(g)}")