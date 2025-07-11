# A python script to display the derived electric field equations.
# The formulas are derived by solving Laplace's equation with the given boundary conditions.

# Define the components for the electric field in Region 1 (0 < phi < pi/2)
# The expression includes the derived number '2' in the numerator and 'pi' in the denominator.
e1_numerator = "2 * sigma_2 * V_0"
e1_denominator = "r * pi * (sigma_1 + sigma_2)"
direction = "i_phi"

# Define the components for the electric field in Region 2 (pi/2 < phi < pi)
e2_numerator = "2 * sigma_1 * V_0"
e2_denominator = "r * pi * (sigma_1 + sigma_2)"

print("The electric field in each region of the resistor is derived as follows:")
print("-" * 60)

# Print the final expression for Region 1
print("Region 1 (0 < phi < pi/2):")
print(f"E_1 = ( {e1_numerator} ) / ( {e1_denominator} ) * {direction}")

print("") # Add a newline for spacing

# Print the final expression for Region 2
print("Region 2 (pi/2 < phi < pi):")
print(f"E_2 = ( {e2_numerator} ) / ( {e2_denominator} ) * {direction}")

print("-" * 60)