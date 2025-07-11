# Parameters for the maximal configuration
# Genus of the fixed curve on the K3 surface S
g = 10
# Number of isolated fixed points on the K3 surface S
k = 18
# Number of fixed points of the involution on the curve C of genus 2
n = 6

# First, calculate the rank of the invariant lattice r_rho
r_rho = 11 - g + k / 2

# Calculate the Hodge number h^{1,1} of the resolved manifold M
# The formula is h¹¹(M) = (r_rho + 1) + n * (1 + k)
h_1_1_M = (r_rho + 1) + n * (1 + k)

# An alternative, more direct formula derived in the thinking steps is h¹¹(M) = 18 - g + 6.5k
# Let's verify with both.
h_1_1_M_direct = 18 - g + 6.5 * k

print("Step 1: Determine the optimal parameters.")
print(f"The involution on the genus 2 curve C is chosen to be hyperelliptic, giving n = {n} fixed points.")
print(f"The involution on the K3 surface S is chosen to maximize k while satisfying constraints, giving g = {g} and k = {k}.")
print("-" * 20)
print("Step 2: Calculate intermediate values.")
print(f"The rank of the invariant lattice, r_rho = 11 - g + k/2 = 11 - {g} + {k}/2 = {int(r_rho)}")
print("-" * 20)
print("Step 3: Calculate the final Hodge number h^{1,1}(M).")
print("Using the formula h¹¹(M) = r_rho + 1 + n * (1 + k):")
print(f"h¹¹(M) = {int(r_rho)} + 1 + {n} * (1 + {k}) = {int(h_1_1_M)}")
print("\nVerifying with the simplified formula h¹¹(M) = 18 - g + 6.5k:")
print(f"h¹¹(M) = 18 - {g} + 6.5 * {k} = {int(h_1_1_M_direct)}")
print("-" * 20)
print(f"The maximal value of the Hodge number h¹¹(M) is: {int(h_1_1_M)}")
