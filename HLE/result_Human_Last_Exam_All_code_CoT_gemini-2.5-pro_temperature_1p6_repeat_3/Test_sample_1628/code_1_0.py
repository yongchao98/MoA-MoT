import math

# --- User-defined Section Properties (in consistent units, e.g., mm) ---
# Web properties
h = 200.0   # Height of the web
t_w = 10.0  # Thickness of the web

# Top flange properties
b_1 = 100.0 # Width of the top flange
t_1 = 12.0  # Thickness of the top flange

# Bottom flange properties
b_2 = 80.0  # Width of the bottom flange
t_2 = 12.0  # Thickness of the bottom flange

# --- Calculations (based on thin-walled section theory) ---

# 1. Calculate Moment of Inertia (I_x) about the horizontal centroidal axis
# This formula is an approximation assuming the centroid is at h/2.
i_x_web = (t_w * h**3) / 12
i_x_flange1 = t_1 * b_1 * (h / 2)**2
i_x_flange2 = t_2 * b_2 * (h / 2)**2
I_x = i_x_web + i_x_flange1 + i_x_flange2

# 2. Calculate the shear center offset 'e' from the web's centerline
# This formula gives the horizontal distance from the web centerline to the shear center.
numerator = h**2 * (t_1 * b_1**2 + t_2 * b_2**2)
denominator = 8 * I_x
e = numerator / denominator

# --- Output the Results ---
print("--- Asymmetric Channel Section Shear Center Calculation ---")
print(f"Given dimensions (mm):")
print(f"  Web Height (h) = {h}, Web Thickness (t_w) = {t_w}")
print(f"  Top Flange (b1, t1) = ({b_1}, {t_1})")
print(f"  Bottom Flange (b2, t2) = ({b_2}, {t_2})\n")

print("--- Step 1: Calculate Moment of Inertia (I_x) ---")
print("Formula (approximate): I_x = (t_w * h^3 / 12) + (t_1 * b_1 * (h/2)^2) + (t_2 * b_2 * (h/2)^2)")
print(f"I_x = ({t_w} * {h}^3 / 12) + ({t_1} * {b_1} * ({h}/2)^2) + ({t_2} * {b_2} * ({h}/2)^2)")
print(f"I_x = {i_x_web:.2f} + {i_x_flange1:.2f} + {i_x_flange2:.2f}")
print(f"Total I_x = {I_x:.2f} mm^4\n")

print("--- Step 2: Calculate Shear Center offset (e) ---")
print("Formula: e = [h^2 * (t_1*b_1^2 + t_2*b_2^2)] / (8 * I_x)")
print(f"e = [{h}^2 * ({t_1}*{b_1}^2 + {t_2}*{b_2}^2)] / (8 * {I_x:.2f})")
print(f"e = [{numerator:.2f}] / [{denominator:.2f}]")
print(f"\nThe shear center is located at a distance e = {e:.2f} mm from the web centerline.")
