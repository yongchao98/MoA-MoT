import math

# Plan:
# 1. Define the geometric properties of an example asymmetric channel section.
# 2. Calculate the position of the section's neutral axis (vertical centroid, y_c).
# 3. Calculate the moment of inertia (I_x) about the neutral axis.
# 4. Use an engineering formula to find the shear center offset 'e' from the web centerline.
# 5. Print all variables and the final result to show the calculation process.

# 1. Define section dimensions (example values in mm)
h = 200.0   # Web height
t_w = 8.0   # Web thickness
b1 = 100.0  # Top flange width
t_f1 = 12.0 # Top flange thickness
b2 = 80.0   # Bottom flange width
t_f2 = 10.0 # Bottom flange thickness

print(f"--- Input Section Dimensions (mm) ---")
print(f"Web height (h): {h}")
print(f"Web thickness (t_w): {t_w}")
print(f"Top flange width (b1): {b1}")
print(f"Top flange thickness (t_f1): {t_f1}")
print(f"Bottom flange width (b2): {b2}")
print(f"Bottom flange thickness (t_f2): {t_f2}\n")


# 2. Calculate the Neutral Axis (y_c)
# Let origin (y=0) be at the outside edge (bottom) of the bottom flange.
A_web = h * t_w
y_web_centroid = t_f2 + h / 2

A_f1 = b1 * t_f1 # Top flange area
y_f1_centroid = t_f2 + h + t_f1 / 2

A_f2 = b2 * t_f2 # Bottom flange area
y_f2_centroid = t_f2 / 2

total_area = A_web + A_f1 + A_f2
# Vertical centroid (y_c) calculation
y_c = (A_web * y_web_centroid + A_f1 * y_f1_centroid + A_f2 * y_f2_centroid) / total_area

print(f"--- Intermediate Calculations ---")
print(f"Total Area: {total_area:.2f} mm^2")
print(f"Vertical centroid (y_c) from bottom edge: {y_c:.2f} mm\n")


# 3. Calculate Moment of Inertia (I_x) using the Parallel Axis Theorem
# I_x = Σ (I_c_i + A_i * d_i^2)
I_web = (t_w * h**3) / 12 + A_web * (y_web_centroid - y_c)**2
I_f1 = (b1 * t_f1**3) / 12 + A_f1 * (y_f1_centroid - y_c)**2
I_f2 = (b2 * t_f2**3) / 12 + A_f2 * (y_f2_centroid - y_c)**2
I_x = I_web + I_f1 + I_f2

print(f"Moment of inertia about neutral axis (I_x): {I_x:.2e} mm^4\n")

# 4. Calculate the shear center offset 'e' from the web's vertical centerline
# The formula used is a common approximation for thin-walled sections.
# e = (h / (4 * I_x)) * (t_f1 * h1 * b1^2 + t_f2 * h2 * b2^2)
# where h1 and h2 are the distances from the neutral axis to the centroids of the flanges.
h1 = y_f1_centroid - y_c
h2 = y_c - y_f2_centroid

numerator = h * (t_f1 * h1 * b1**2 + t_f2 * h2 * b2**2)
denominator = 4 * I_x
e = numerator / denominator

# 5. Print the breakdown of the final calculation
print(f"--- Shear Center Calculation ---")
print(f"This calculation finds the horizontal offset 'e' of the shear center from the web's centerline.")
print(f"Using formula: e = (h / (4 * I_x)) * (t_f1 * h1 * b1^2 + t_f2 * h2 * b2^2)")
print("\nValues used in the equation:")
print(f"h = {h:.2f}")
print(f"I_x = {I_x:.2e}")
print(f"t_f1 = {t_f1:.2f}, h1 = {h1:.2f}, b1 = {b1:.2f}")
print(f"t_f2 = {t_f2:.2f}, h2 = {h2:.2f}, b2 = {b2:.2f}")

print("\nFinal equation with numbers:")
print(f"e = ({h:.2f} / (4 * {I_x:.2e})) * ({t_f1:.2f} * {h1:.2f} * {b1**2:.2f} + {t_f2:.2f} * {h2:.2f} * {b2**2:.2f})")
print(f"e = {e:.2f} mm")
