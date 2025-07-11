import math

# This script calculates the shear center for an asymmetric channel section.
# Plan:
# 1. Define the dimensions of an example asymmetric channel section.
# 2. Calculate the location of the horizontal centroidal axis (Neutral Axis for bending about x-axis).
# 3. Calculate the moment of inertia (I_x) about this axis using the parallel axis theorem.
# 4. Use the derived formula for shear center eccentricity 'e' based on shear flow principles.
# 5. Print all results and calculation steps clearly.

# --- Step 1: Define the dimensions of the asymmetric channel section ---
# (All dimensions in mm)
print("--- Input Section Dimensions ---")
b1 = 50.0   # Top flange width
t1 = 10.0   # Top flange thickness
b2 = 100.0  # Bottom flange width
t2 = 12.0   # Bottom flange thickness
h_web = 180.0 # Clear height of the web
tw = 8.0    # Web thickness
print(f"Top Flange (b1 x t1): {b1} mm x {t1} mm")
print(f"Bottom Flange (b2 x t2): {b2} mm x {t2} mm")
print(f"Web (h_web x tw): {h_web} mm x {tw} mm")


# --- Step 2: Calculate Centroid (y_bar) ---
# Calculate areas of the three components
A1 = b1 * t1      # Top flange
A2 = h_web * tw   # Web
A3 = b2 * t2      # Bottom flange
A_total = A1 + A2 + A3

# Calculate the vertical centroid (y_bar) from the bottom outer edge
# y-coordinates of the centroid of each component from the bottom edge
y1_c = h_web + t2 + t1 / 2.0 # Centroid of top flange
y2_c = t2 + h_web / 2.0      # Centroid of web
y3_c = t2 / 2.0              # Centroid of bottom flange

# The final equation for the centroid y_bar
y_bar_numerator = (A1 * y1_c) + (A2 * y2_c) + (A3 * y3_c)
y_bar = y_bar_numerator / A_total

print("\n--- Centroid Calculation (from bottom edge) ---")
print(f"y_bar = (A_top*y_top + A_web*y_web + A_bottom*y_bottom) / (A_total)")
print(f"y_bar = (({A1:.1f}*{y1_c:.1f}) + ({A2:.1f}*{y2_c:.1f}) + ({A3:.1f}*{y3_c:.1f})) / ({A_total:.1f})")
print(f"y_bar = {y_bar_numerator:.1f} / {A_total:.1f} = {y_bar:.2f} mm")


# --- Step 3: Calculate the Moment of Inertia (I_x) ---
# Using the Parallel Axis Theorem: I_x = sum(I_xc + A * d^2)
# I_xc for each rectangular component
I_xc1 = (b1 * t1**3) / 12.0
I_xc2 = (tw * h_web**3) / 12.0
I_xc3 = (b2 * t2**3) / 12.0

# Distances (d) from component centroids to the section's centroidal axis
d1 = y1_c - y_bar      # Distance for top flange
d2 = y2_c - y_bar      # Distance for web
d3 = y_bar - y3_c      # Positive distance for bottom flange

# Individual terms for the parallel axis theorem
I_x1 = I_xc1 + A1 * d1**2
I_x2 = I_xc2 + A2 * d2**2
I_x3 = I_xc3 + A3 * d3**2
I_x = I_x1 + I_x2 + I_x3

print("\n--- Moment of Inertia Calculation (I_x) ---")
print("I_x = sum(I_xc + A*d^2)")
print(f"I_x_top_flange    = {I_xc1:.1f} + {A1:.1f}*({d1:.2f})^2 = {I_x1:.2f} mm^4")
print(f"I_x_web           = {I_xc2:.1f} + {A2:.1f}*({d2:.2f})^2 = {I_x2:.2f} mm^4")
print(f"I_x_bottom_flange = {I_xc3:.1f} + {A3:.1f}*({d3:.2f})^2 = {I_x3:.2f} mm^4")
print(f"I_x_total = {I_x1:.2f} + {I_x2:.2f} + {I_x3:.2f} = {I_x:.2f} mm^4")


# --- Step 4: Calculate the Shear Center eccentricity (e) ---
# The formula for eccentricity 'e' from the web's centerline is:
# e = (1 / (2*I_x)) * [ (b1^2 * t1 * d1) + (b2^2 * t2 * d2) ]
# where d1 and d2 are the distances from the neutral axis to the
# centerlines of the top and bottom flanges, respectively.

e_numerator = (b1**2 * t1 * d1) + (b2**2 * t2 * d3)
e_denominator = 2 * I_x
e = e_numerator / e_denominator

print("\n--- Shear Center Calculation ---")
print("e = [ (b1^2*t1*d_top) + (b2^2*t2*d_bottom) ] / (2 * I_x)")
print(f"Numerator = ({b1**2:.1f}*{t1:.1f}*{d1:.2f}) + ({b2**2:.1f}*{t2:.1f}*{d3:.2f}) = {e_numerator:.2f}")
print(f"Denominator = 2 * {I_x:.2f} = {e_denominator:.2f}")
print(f"e = {e_numerator:.2f} / {e_denominator:.2f} = {e:.2f} mm")

print("\n--- Final Result ---")
print("The shear center is the point where a shear force can be applied without causing twisting.")
print(f"For this asymmetric channel, the shear center is located {e:.2f} mm horizontally from the web's centerline.")
print("This location is outside the section, on the side opposite the flanges.")
