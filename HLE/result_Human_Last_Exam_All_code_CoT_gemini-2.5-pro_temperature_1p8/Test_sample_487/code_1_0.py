# Parameters for the involution on the curve C of genus 2
# We choose the hyperelliptic involution to maximize the number of fixed points.
n = 6  # Number of fixed points for psi on C

# Parameters for the non-symplectic involution on the K3 surface S
# To maximize h^{1,1}(M) = 21 - r + n*k, we analyze the possible (r, k) pairs.
# The maximum is found for the family where the resolved quotient is a K3 surface.
# In this family, h^{1,1}(M) = 93 - 7*r. To maximize, we need to minimize r.
# The minimum possible rank for the anti-invariant lattice in this family is r=2.
r = 2  # Rank of the anti-invariant lattice for rho on S
k = 12 - r  # Number of fixed rational curves for rho on S

# Calculate the Hodge number of the quotient orbifold Y = (S x C) / sigma
h11_Y = 21 - r

# Calculate the contribution from blowing up the n*k singular curves
blowup_contribution = n * k

# Calculate the final Hodge number for the resolved manifold M
h11_M = h11_Y + blowup_contribution

# Print the final result, showing each number in the equation
print(f"The maximal value of h^(1,1)(M) is obtained with n={n}, r={r}, k={k}.")
print(f"The equation for h^(1,1)(M) is: h^(1,1) = (21 - r) + (n * k)")
print(f"Plugging in the values: {h11_M} = (21 - {r}) + ({n} * {k})")

print(f"Final Answer is {h11_M}")
<<<79>>>