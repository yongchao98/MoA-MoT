# Define the known variables based on the problem description.
# The problem uses 'm^2' as a symbolic value for the squared mass
# of the 5 tensor degrees of freedom. We represent it as a string.
m_tensor_squared_str = "m^2"
num_tensor_dof = 5
num_scalar_dof = 1

# From the analysis of the equations of motion for this specific theory,
# a key relationship is found between the squared mass of the tensor modes (M_T^2)
# and the squared mass of the scalar mode (M_S^2).
# The relation is M_T^2 = 2 * M_S^2.
ratio_tensor_to_scalar_mass_sq = 2

# We are given the squared mass of the 5 tensor modes.
print(f"The squared mass of the {num_tensor_dof} tensor modes is given as: M_T^2 = {m_tensor_squared_str}")

# Print the relationship between the masses.
print(f"The derived relationship between the squared masses is: M_T^2 = {ratio_tensor_to_scalar_mass_sq} * M_S^2")

# Now, we solve for the squared mass of the 6th degree of freedom (the scalar mode).
# M_S^2 = M_T^2 / 2
m_scalar_squared_str = f"{m_tensor_squared_str} / {ratio_tensor_to_scalar_mass_sq}"

# Print the final result for the squared mass of the sixth degree of freedom.
print(f"\nBy substituting the given information, we can calculate the squared mass of the {num_scalar_dof+num_tensor_dof}th degree of freedom:")
print(f"M_S^2 = {m_tensor_squared_str} / {ratio_tensor_to_scalar_mass_sq}")
print(f"The squared mass of the sixth degree of freedom is {m_scalar_squared_str}.")
