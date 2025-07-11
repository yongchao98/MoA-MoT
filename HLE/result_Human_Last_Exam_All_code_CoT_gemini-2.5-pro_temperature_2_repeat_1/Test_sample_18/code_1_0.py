import numpy as np

# Step 1: Define the K-matrix for the bosons.
# The problem states that the K-matrix for the nu=2 Bosonic Integer Quantum Hall state
# is the Pauli matrix sigma_x.
K_boson = np.array([[0, 1], 
                    [1, 0]])

# Let's verify the filling fraction. For a K-matrix, nu = q^T * K^-1 * q.
# For sigma_x, its inverse is itself: K_boson_inv = K_boson.
# Assuming unit charge bosons, the charge vector q is [1, 1].
# nu = [1, 1] * [[0, 1], [1, 0]] * [[1], [1]] = [1, 1] * [[1], [1]] = 2.
# The given K_matrix is consistent with nu=2.

# Step 2: Determine the K-matrix for the constituent composite fermions.
# The bosons are Cooper pairs of composite fermions (b ~ cf^2). This relationship implies
# the K-matrix for the composite fermions (K_CF) is 4 times the K-matrix of the bosons.
K_cf_raw = 4 * K_boson

# Step 3: Enforce the fermion statistics for the composite fermions.
# K_cf_raw is [[0, 4], [4, 0]]. Its diagonal elements are 0 (even), which corresponds to bosons.
# However, the problem specifies these are Composite *Fermions*. A K-matrix for fermions
# must have odd integer diagonal elements. The simplest way to satisfy this is to assume
# the composite fermions form a base state of nu=1 for each species, which corresponds to
# adding the identity matrix I to K_cf_raw.
identity_matrix = np.identity(2, dtype=int)
K_cf = K_cf_raw + identity_matrix

# Step 4: Undo the flux attachment to find the K-matrix of the original fermions.
# The composite fermions are described as having two flux quanta attached (m=2).
# Attaching 'm' fluxes per particle corresponds to adding 2*m to the diagonal elements
# of the K-matrix. So, K_CF = K_fermion + (2 * m) * I.
# With m=2, this is K_CF = K_fermion + 4 * I.
# To find the K-matrix of the original fermions, K_fermion, we rearrange the formula:
# K_fermion = K_CF - 4 * I.
m = 2
K_fermion = K_cf - (2 * m) * identity_matrix

# Step 5: Print the final result.
# The final K-matrix describes the fractional state of the underlying fermions.
print("The plan is as follows:")
print("1. Start with the given K_boson = sigma_x.")
print(f"K_boson = \n{K_boson}\n")
print("2. The bosons are Cooper-pairs of composite fermions, so K_CF is related to 4 * K_boson.")
print("3. We enforce fermion statistics on K_CF by ensuring its diagonals are odd.")
print(f"K_CF = 4 * K_boson + I = \n{K_cf}\n")
print("4. We find the original fermion K-matrix, K_f, by reversing the attachment of m=2 fluxes.")
print(f"K_f = K_CF - 2*m*I = K_CF - 4*I\n")
print("The final K-matrix for the resulting fractional state is:")

# We output each number in the final equation as requested.
print(f"K = [ [ {K_fermion[0, 0]}, {K_fermion[0, 1]} ],")
print(f"      [ {K_fermion[1, 0]}, {K_fermion[1, 1]} ] ]")
