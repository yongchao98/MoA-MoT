# Part A: Fiber Dimension Calculation
# The physical field is the torsion tensor T^k_{ij}, which is required for a metric-compatible connection.
# The torsion tensor is antisymmetric in its lower two indices (i, j).
# For each upper index k (from 1 to 3), T^k is a 3x3 antisymmetric matrix.
# A 3x3 antisymmetric matrix has 3 independent components.
num_upper_indices = 3
components_per_antisymmetric_matrix = 3
dimension_of_fiber = num_upper_indices * components_per_antisymmetric_matrix

# Part B: Energy Coefficients Calculation
# The energy E is a quadratic function of the torsion tensor components.
# The number of coefficients equals the number of independent quadratic invariants of a
# general second-rank tensor under cubic symmetry.
# We find these by decomposing the tensor into its symmetric (S) and antisymmetric (A) parts.

# For a symmetric 3x3 tensor under cubic symmetry, there are 3 quadratic invariants.
num_invariants_symmetric = 3
# For an antisymmetric 3x3 tensor (an axial vector) under cubic symmetry, there is 1 quadratic invariant.
num_invariants_antisymmetric = 1

# The total number of coefficients is the sum of the invariants from the symmetric and antisymmetric parts.
num_energy_coefficients = num_invariants_symmetric + num_invariants_antisymmetric

# Print the final answer in the requested format "A B"
print(f"{dimension_of_fiber} {num_energy_coefficients}")
