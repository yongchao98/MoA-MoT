import math

# Step 1: Define the geometry of an asymmetric channel section (in mm)
# An asymmetric channel with a taller top flange
h = 150.0   # Height of the web
b1 = 100.0  # Width of the top flange
tf1 = 12.0  # Thickness of the top flange
b2 = 75.0   # Width of the bottom flange
tf2 = 10.0  # Thickness of the bottom flange
tw = 8.0    # Thickness of the web

print("--- Asymmetric Channel Section Properties (mm) ---")
print(f"Web height (h): {h}")
print(f"Top flange width (b1): {b1}")
print(f"Top flange thickness (tf1): {tf1}")
print(f"Bottom flange width (b2): {b2}")
print(f"Bottom flange thickness (tf2): {tf2}")
print(f"Web thickness (tw): {tw}\n")


# Step 2: Calculate the vertical centroid (y_bar) from a datum
# Datum is the bottom outer edge of the bottom flange
# The section is composed of three rectangular areas
A1 = b1 * tf1  # Area of top flange
y1_local = tf2 + h + tf1 / 2.0  # Centroid of top flange from datum

A2 = h * tw    # Area of web
y2_local = tf2 + h / 2.0  # Centroid of web from datum

A3 = b2 * tf2  # Area of bottom flange
y3_local = tf2 / 2.0  # Centroid of bottom flange from datum

total_area = A1 + A2 + A3
sum_Ay = (A1 * y1_local) + (A2 * y2_local) + (A3 * y3_local)
y_bar = sum_Ay / total_area

print("--- Centroid Calculation ---")
print(f"Total Area = {total_area:.2f} mm^2")
print(f"Vertical Centroid (y_bar) from bottom edge = {y_bar:.2f} mm\n")


# Step 3: Calculate the Moment of Inertia (I_x) about the horizontal centroidal axis
# Using the Parallel Axis Theorem: I = I_c + A*d^2
# For top flange
Ic1 = (b1 * tf1**3) / 12.0
d1 = y1_local - y_bar
Ix1 = Ic1 + A1 * d1**2

# For web
Ic2 = (tw * h**3) / 12.0
d2 = y2_local - y_bar
Ix2 = Ic2 + A2 * d2**2

# For bottom flange
Ic3 = (b2 * tf2**3) / 12.0
d3 = y3_local - y_bar
Ix3 = Ic3 + A3 * d3**2

I_x_total = Ix1 + Ix2 + Ix3

print("--- Moment of Inertia Calculation ---")
print(f"Moment of Inertia (I_x) = {I_x_total:.2f} mm^4\n")


# Step 4: Calculate the Shear Center location (e)
# 'e' is the horizontal distance from the web centerline.
# The formula is e = (1 / (2*I_x)) * [tf1*h1^2*b1^2 + tf2*h2^2*b2^2]
# where h1 and h2 are the distances from the centroid to the centerlines of the flanges.
h1 = d1  # This is the same as the distance 'd' for the top flange
h2 = abs(d3) # This is the same as the distance 'd' for the bottom flange

term1 = tf1 * h1**2 * b1**2
term2 = tf2 * h2**2 * b2**2
e_numerator = term1 + term2
e = e_numerator / (2 * I_x_total)

print("--- Shear Center Calculation ---")
print("Formula: e = [tf1*h1^2*b1^2 + tf2*h2^2*b2^2] / (2 * I_x)")
print("where h1 and h2 are the distances from the section centroid to the flange centerlines.")
print(f"h1 = |y1_local - y_bar| = |{y1_local:.2f} - {y_bar:.2f}| = {h1:.2f} mm")
print(f"h2 = |y_bar - y3_local| = |{y_bar:.2f} - {y3_local:.2f}| = {h2:.2f} mm")
print("\nFinal Equation with values:")
# Final printout as requested, showing the full equation with numbers
print(f"e = [{tf1:.1f}*{h1:.2f}^2*{b1:.1f}^2 + {tf2:.1f}*{h2:.2f}^2*{b2:.1f}^2] / (2 * {I_x_total:.2f})")
print(f"e = [{term1:.2f} + {term2:.2f}] / (2 * {I_x_total:.2f})")
print(f"e = {e_numerator:.2f} / {2*I_x_total:.2f}")

print(f"\nResult: The shear center is located at a distance e = {e:.2f} mm from the web's centerline.")
print("This location is outside the physical cross-section.")
