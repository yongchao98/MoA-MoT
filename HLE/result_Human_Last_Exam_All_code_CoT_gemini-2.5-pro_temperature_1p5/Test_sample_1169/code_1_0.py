# The number of initial data required is the dimension of the physical phase space.
# We found this by counting the independent initial coordinates and momenta.
#
# Number of independent initial coordinates:
# x_mu has 4 components.
# w_mu has 4 components with one constraint (w^2=1), so 3 independent components.
# Total coordinate data = 4 + 3 = 7
num_coord_data = 7

# Number of independent initial momenta:
# p_x_mu has 4 components.
# p_w_mu has 4 components.
# Total momenta components = 4 + 4 = 8
#
# Constraints on momenta:
# 1. Hamiltonian constraint H=0 from reparametrization invariance.
# 2. p_x . w = 0 from the form of the Lagrangian.
# 3. p_w . w = 0 from the form of the Lagrangian.
# Total constraints = 3
# Total momentum data = 8 - 3 = 5
num_momentum_data = 5

# Total number of initial data points
total_initial_data = num_coord_data + num_momentum_data

print(total_initial_data)