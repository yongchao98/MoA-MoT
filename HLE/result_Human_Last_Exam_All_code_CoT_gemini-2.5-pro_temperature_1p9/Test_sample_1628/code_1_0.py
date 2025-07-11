import math

# Step 1: Define the geometry of an asymmetric channel section
# All units are in millimeters (mm)
h = 150.0   # Web height
t_w = 10.0  # Web thickness
b_1 = 100.0 # Top flange width
t_1 = 15.0  # Top flange thickness
b_2 = 75.0  # Bottom flange width
t_2 = 12.0  # Bottom flange thickness

print("--- Geometric Properties ---")
print(f"Web height (h): {h} mm")
print(f"Web thickness (t_w): {t_w} mm")
print(f"Top flange width (b_1): {b_1} mm")
print(f"Top flange thickness (t_1): {t_1} mm")
print(f"Bottom flange width (b_2): {b_2} mm")
print(f"Bottom flange thickness (t_2): {t_2} mm\n")

# Step 2: Calculate the vertical centroid (c_y) from the bottom of the section
# Areas of the components
A_1 = b_1 * t_1  # Top flange
A_2 = b_2 * t_2  # Bottom flange
A_web = h * t_w  # Web
A_total = A_1 + A_2 + A_web

# y-coordinates of the centroid of each part from the bottom edge
y_1 = t_2 + h + t_1 / 2.0
y_2 = t_2 / 2.0
y_web = t_2 + h / 2.0

# Calculate vertical centroid c_y
c_y = (A_1 * y_1 + A_2 * y_2 + A_web * y_web) / A_total

print("--- Centroid Calculation ---")
print(f"Vertical centroid (c_y) from bottom edge: {c_y:.2f} mm\n")

# Step 3: Calculate Moment of Inertia (I_x) about the horizontal centroidal axis
# Using the Parallel Axis Theorem: I_x = Σ(I_c + A*d^2)
# I_c for a rectangle is (base * height^3) / 12

# Top flange contribution to I_x
I_x_1 = (b_1 * t_1**3) / 12.0 + A_1 * (y_1 - c_y)**2
# Bottom flange contribution to I_x
I_x_2 = (b_2 * t_2**3) / 12.0 + A_2 * (y_2 - c_y)**2
# Web contribution to I_x
I_x_web = (t_w * h**3) / 12.0 + A_web * (y_web - c_y)**2

I_x_total = I_x_1 + I_x_2 + I_x_web

print("--- Moment of Inertia Calculation ---")
print(f"Moment of inertia (I_x) about centroidal axis: {I_x_total:.2f} mm^4\n")

# Step 4: Calculate the shear center offset (e) from the web centerline
# Distances from the neutral axis (at c_y) to the flange centerlines
d_1 = y_1 - c_y
d_2 = c_y - y_2

# Numerator and denominator of the shear center formula
numerator = (d_1 * t_1 * b_1**2 + d_2 * t_2 * b_2**2)
denominator = 2 * I_x_total
e = numerator / denominator

print("--- Shear Center Calculation ---")
print("The formula for the shear center offset 'e' from the web centerline is:")
print("e = (d_1 * t_1 * b_1^2 + d_2 * t_2 * b_2^2) / (2 * I_x)\n")

print("Plugging in the calculated values:")
print(f"d_1 (distance from NA to top flange centerline) = {d_1:.2f} mm")
print(f"d_2 (distance from NA to bottom flange centerline) = {d_2:.2f} mm\n")

# As per the user request, printing the equation with all the numbers.
print("Final Equation:")
print(f"e = (({d_1:.2f} * {t_1:.2f} * {b_1:.2f}^2) + ({d_2:.2f} * {t_2:.2f} * {b_2:.2f}^2)) / (2 * {I_x_total:.2f})")

# Calculate the value of each term in the numerator for clarity
term1_numerator = d_1 * t_1 * b_1**2
term2_numerator = d_2 * t_2 * b_2**2
print(f"e = ({term1_numerator:.2f} + {term2_numerator:.2f}) / {denominator:.2f}")
print(f"e = {numerator:.2f} / {denominator:.2f}\n")

print("--- Result ---")
print(f"The shear center is located at a distance 'e' = {e:.2f} mm horizontally from the web centerline.")
print("This position is outside the physical cross-section, on the side away from the flanges.")