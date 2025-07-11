import math

# This script calculates the shear center location for an asymmetric channel section.

# 1. Define the geometry of the asymmetric channel section (dimensions in mm)
b1 = 100.0  # Top flange width
t_f1 = 15.0 # Top flange thickness
b2 = 75.0   # Bottom flange width
t_f2 = 12.0 # Bottom flange thickness
h_w = 200.0 # Web height (clear distance between flanges)
t_w = 10.0  # Web thickness

print("Calculating shear center for an asymmetric channel section with the following dimensions (mm):")
print(f"Top Flange: width(b1)={b1}, thickness(t_f1)={t_f1}")
print(f"Bottom Flange: width(b2)={b2}, thickness(t_f2)={t_f2}")
print(f"Web: height(h_w)={h_w}, thickness(t_w)={t_w}\n")

# 2. Calculate the area of each component
A1 = b1 * t_f1   # Area of top flange
A_w = h_w * t_w  # Area of web
A2 = b2 * t_f2   # Area of bottom flange
A_total = A1 + A_w + A2

# 3. Find the vertical centroid (y_bar) of the section.
# We establish a reference axis at the bottom edge of the bottom flange.
# First, find the y-coordinate of the centroid of each part.
y1_c = t_f2 + h_w + t_f1 / 2.0  # Centroid of top flange
y_w_c = t_f2 + h_w / 2.0        # Centroid of web
y2_c = t_f2 / 2.0               # Centroid of bottom flange

# Calculate the weighted average to find the section's centroidal axis location.
y_bar_numerator = (A1 * y1_c + A_w * y_w_c + A2 * y2_c)
y_bar = y_bar_numerator / A_total
print(f"The vertical centroid (y_bar) from the bottom edge is: {y_bar:.2f} mm\n")

# 4. Calculate the Moment of Inertia (I_x) about the centroidal x-axis using the Parallel Axis Theorem.
# I_x = Σ (I_c + A*d^2) for each part
I_x1 = (b1 * t_f1**3) / 12.0 + A1 * (y1_c - y_bar)**2
I_xw = (t_w * h_w**3) / 12.0 + A_w * (y_w_c - y_bar)**2
I_x2 = (b2 * t_f2**3) / 12.0 + A2 * (y2_c - y_bar)**2
I_x = I_x1 + I_xw + I_x2
print(f"The moment of inertia (I_x) about the centroidal axis is: {I_x:.2f} mm^4\n")

# 5. Calculate the distances d1 and d2 from the section centroid (y_bar) to the centroids of the top and bottom flanges.
d1 = y1_c - y_bar
d2 = abs(y2_c - y_bar)
print(f"Distance from Neutral Axis to Top Flange centroid (d1): {d1:.2f} mm")
print(f"Distance from Neutral Axis to Bottom Flange centroid (d2): {d2:.2f} mm\n")

# 6. Calculate the shear center offset 'e_c' from the vertical centroidal axis.
# The formula is derived from balancing the external torque (V*e) with the internal torque from flange shear flows.
# e_c = (1 / (2 * I_x)) * (t_f1 * b1^2 * d1^2 + t_f2 * b2^2 * d2^2)
numerator_e_c = (t_f1 * b1**2 * d1**2 + t_f2 * b2**2 * d2**2)
e_c = numerator_e_c / (2 * I_x)

print("The formula for the shear center offset from the centroid (e_c) is:")
print("e_c = (t_f1*b1^2*d1^2 + t_f2*b2^2*d2^2) / (2 * I_x)\n")

print("Final Calculation:")
# Print the final calculation with all the numbers plugged in.
print(f"e_c = ({t_f1:.2f}*{b1:.2f}²*{d1:.2f}² + {t_f2:.2f}*{b2:.2f}²*{d2:.2f}²) / (2 * {I_x:.2f})")
print(f"e_c = {numerator_e_c:.2f} / {2*I_x:.2f}")
print(f"e_c = {e_c:.2f} mm")
print("\nThe shear center is located {:.2f} mm from the vertical centroidal axis, on the side of the web (away from the flanges).".format(e_c))