import math

# Step 1: Define the problem based on the physical model.
# The string consists of two parts:
# 1. A quarter-circular arc of radius R.
# 2. A vertical segment hanging from the end of the arc.
# The length of the hanging part (L_h) is determined by tension balance.
# Tension at the end of the arc (theta=pi/2) is T = rho * g * R * (1 - cos(pi/2)) = rho * g * R.
# Weight of hanging part is W_h = (rho * L_h) * g.
# For equilibrium, T = W_h, which means rho * g * R = rho * L_h * g => L_h = R.

# Step 2: Define properties of each part.
# Let R=1 for simplicity to find the numerical coefficients.
# Let rho=1 for simplicity as it cancels out.
R = 1.0

# Part 1: Quarter-circular arc (from theta=0 to pi/2)
# Mass M1 = rho * length = 1 * (pi * R / 2) = pi/2
M1 = math.pi / 2
# Center of mass of a quarter-circular arc is (2R/pi, 2R/pi).
X1 = 2 * R / math.pi
Z1 = 2 * R / math.pi

# Part 2: Vertical hanging segment
# It hangs from (R, 0) down to (R, -R).
# Mass M2 = rho * length = 1 * R = 1
M2 = R
# Center of mass of the vertical segment is its midpoint.
X2 = R
Z2 = -R / 2

# Step 3: Calculate the overall center of mass (X_cm, Z_cm).
# M_total = M1 + M2
# X_cm = (M1*X1 + M2*X2) / M_total
# Z_cm = (M1*Z1 + M2*Z2) / M_total
M_total = M1 + M2

X_cm_numerator = M1 * X1 + M2 * X2
# Analytically: (pi/2)*(2R/pi) + R*R = R + R = 2R. With R=1, this is 2.
Z_cm_numerator = M1 * Z1 + M2 * Z2
# Analytically: (pi/2)*(2R/pi) + R*(-R/2) = R - R/2 = R/2. With R=1, this is 0.5.

# The final coordinates are X_cm = (2R) / (pi/2 + 1) = 4R / (pi + 2)
# and Z_cm = (R/2) / (pi/2 + 1) = R / (pi + 2).
# The requested "raw numbers" are the coefficients of R.
k_x_numerator = 4
k_z_numerator = 1
denominator = math.pi + 2

# Step 4: Calculate the numerical coefficients and print results.
k_x = k_x_numerator / denominator
k_z = k_z_numerator / denominator

print(f"The equation for the horizontal coordinate coefficient is: {k_x_numerator} / (pi + 2)")
print(f"The equation for the vertical coordinate coefficient is: {k_z_numerator} / (pi + 2)")
print("The raw numbers for the horizontal and vertical coordinates (coefficients of R) are:")
print(f"{k_x},{k_z}")