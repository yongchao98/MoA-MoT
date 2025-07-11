import math

# Step 1: Define the problem parameters based on the physical analysis.
# The string consists of two parts in equilibrium:
# 1. An arc on the pumpkin: A quarter-circle of radius R.
# 2. A hanging part: A vertical line of length R.
# The linear mass density is rho. The radius of the pumpkin is R.

# Step 2: Define the properties of each part.
# The mass of a part is rho * length.
# Arc properties:
# Length_arc = (math.pi * R) / 2
# Mass_arc = rho * (math.pi * R) / 2
# CoM_arc = (2*R/math.pi, 2*R/math.pi) for (x, z) coordinates.

# Hanging line properties:
# Length_hang = R
# Mass_hang = rho * R
# CoM_hang = (R, -R/2) for (x, z) coordinates.

# Step 3: Calculate the total center of mass.
# We can cancel out R and rho from the weighted average calculation to find the coefficients.
# x_cm = (Mass_arc * x_cm_arc + Mass_hang * x_cm_hang) / (Mass_arc + Mass_hang)
# x_cm = ( (rho*pi*R/2)*(2*R/pi) + (rho*R)*(R) ) / ( rho*pi*R/2 + rho*R )
# x_cm = ( R**2 + R**2 ) / ( pi*R/2 + R )
# x_cm = 2 * R**2 / (R * (pi/2 + 1))
# x_cm = 2 * R / (pi/2 + 1)
# x_cm = 4 * R / (pi + 2)
# The coefficient for x is 4 / (pi + 2).

# z_cm = (Mass_arc * z_cm_arc + Mass_hang * z_cm_hang) / (Mass_arc + Mass_hang)
# z_cm = ( (rho*pi*R/2)*(2*R/pi) + (rho*R)*(-R/2) ) / ( rho*pi*R/2 + rho*R )
# z_cm = ( R**2 - R**2/2 ) / ( pi*R/2 + R )
# z_cm = (R**2/2) / (R * (pi/2 + 1))
# z_cm = (R/2) / (pi/2 + 1)
# z_cm = R / (pi + 2)
# The coefficient for z is 1 / (pi + 2).

# The final equations for the coordinates are:
# x_cm = (4 / (math.pi + 2)) * R
# z_cm = (1 / (math.pi + 2)) * R

# We will now calculate and print these numerical coefficients.
# The horizontal coordinate coefficient:
x_coeff = 4 / (math.pi + 2)

# The vertical coordinate coefficient:
z_coeff = 1 / (math.pi + 2)

# The prompt asks for the raw numbers of the horizontal and vertical coordinates, separated by a comma.
# It also mentions to "output each number in the final equation". We interpret this as providing
# the calculated numerical coefficients for the x and z coordinates.
print(f"{x_coeff},{z_coeff}")