# The Lindhard polarization function, Pi_0(k, omega), describes the linear response of
# the electron density to an external potential. We are interested in its value
# in the static (omega = 0) and long-wavelength (k -> 0) limit.
#
# In this limit, the theory shows that the Lindhard function is equal to the
# negative of the density of states at the Fermi energy, g(epsilon_F).
# Formula: Pi_0(k=0, omega=0) = -g(epsilon_F)
#
# The density of states g(epsilon_F) is not a universal constant; it depends on the
# electron density of the material. To find a single numerical value, we calculate
# the dimensionless Lindhard function, which is normalized by g(epsilon_F).
# Dimensionless Formula: Pi_tilde_0 = Pi_0 / g(epsilon_F)

# For the purpose of the calculation, we can represent g(epsilon_F) as a variable.
# Since we are calculating a ratio, its actual value does not matter, so we can set it to 1.
density_of_states_g_eF = 1.0

# According to the formula, the Lindhard function value is the negative of g(epsilon_F).
lindhard_function_Pi_0 = -density_of_states_g_eF

# The universal dimensionless value is the ratio of Pi_0 to g(epsilon_F).
dimensionless_value = lindhard_function_Pi_0 / density_of_states_g_eF

# Print out the logic and the numbers involved in the final calculation.
print("The value of the Lindhard function in the static, long-wavelength limit is Pi_0 = -g(epsilon_F).")
print("To get a universal number, we compute the dimensionless value: Pi_tilde_0 = Pi_0 / g(epsilon_F).")
print("\nLet's calculate this using representative numbers:")
print(f"Value representing g(epsilon_F): {density_of_states_g_eF}")
print(f"Resulting value for Pi_0 = -g(epsilon_F): {lindhard_function_Pi_0}")
print("\nFinal equation for the dimensionless value:")
print(f"{lindhard_function_Pi_0} / {density_of_states_g_eF} = {dimensionless_value}")
