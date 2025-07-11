import numpy as np

# We found that the angle alpha is given by arctan(3*sqrt(3)).
# The numbers in the final equation are 3 and 3 (from sqrt(3)).
# We will print the derivation of this exact value based on our calculation.

# Calculation step:
# phi_max is the angle of the complex number z = -2.5 + 7.5j * sqrt(3)
# In our derivation, z_numerator = 12*z = -5/2 + i*15*sqrt(3)/2, so this doesn't affect the angle.
# Let's use the numerator values.
re_part = -5.0 / 2.0
im_part = 15.0 * np.sqrt(3) / 2.0

# Calculate phi_max = atan2(Im, Re)
phi_max_rad = np.arctan2(im_part, re_part)

# Calculate alpha = pi - phi_max
alpha_rad = np.pi - phi_max_rad

# For the output, we express alpha in terms of arctan.
# alpha = arctan(3 * sqrt(3))
# Let's print the numbers in the final equation as requested.
number1 = 3
number2 = 3

print(f"The exact value for the angle alpha is derived from the expression arctan(x), where x is computed from the stability function at theta=pi/3.")
print(f"We calculated the corresponding complex value z's components to be proportional to: Real = {re_part}, Imaginary = {im_part}.")
print(f"The tangent of the angle of z is Im/Re = ({im_part}) / ({re_part}) = {im_part/re_part:.4f}.")
print(f"This leads to phi_max = pi + arctan(-3*sqrt(3)) = pi - arctan(3*sqrt(3)).")
print(f"The stability angle is alpha = pi - phi_max.")
print(f"Therefore, alpha = arctan({number1} * sqrt({number2})).")

# Verification
computed_alpha = np.arctan(3 * np.sqrt(3))
print(f"Calculated alpha (radians): {alpha_rad:.6f}")
print(f"Formula alpha (radians):   {computed_alpha:.6f}")