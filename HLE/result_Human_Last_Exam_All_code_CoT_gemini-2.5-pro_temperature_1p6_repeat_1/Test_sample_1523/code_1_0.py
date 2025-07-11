# In lattice field theory, physical calculations are performed on a discretized
# spacetime lattice to regulate the infinite degrees of freedom of the continuum.
# The expectation value of the inner product (ϕ, D_ϕ) for a free scalar field
# on a lattice is equal to the total number of degrees of freedom, i.e.,
# the number of lattice sites.

# This script calculates this value for a given lattice size.

# Define the dimensions of the 4D spacetime lattice
Nt = 24  # Number of sites in the time direction
Ns = 16  # Number of sites in each of the 3 spatial directions

# The total number of degrees of freedom (N_dof) is the total number of sites.
# N_dof = Nt * Ns * Ns * Ns
N_dof = Nt * (Ns**3)

# The result is expressed as an equation showing the calculation.
# We print the value for each term in the final equation as requested.
print(f"The expectation value of the inner product on a {Nt} x {Ns}^3 lattice is:")
print(f"⟨(ϕ, D_ϕ)⟩ = {Nt} * {Ns}**3 = {N_dof}")
