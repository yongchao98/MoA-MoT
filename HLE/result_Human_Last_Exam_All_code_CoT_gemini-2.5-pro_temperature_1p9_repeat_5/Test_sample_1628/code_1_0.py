import math

# Description: This script calculates the shear center location for a thin-walled,
# asymmetric channel section subjected to a vertical shear force.

# 1. Define the geometric properties of the section (all units in mm)
b1 = 100.0  # Top flange width
b2 = 50.0   # Bottom flange width
h = 200.0   # Height between flange centerlines (web height)
t_f = 10.0  # Flange thickness (assumed constant)
t_w = 8.0   # Web thickness

print("Calculating Shear Center for an Asymmetric Channel Section")
print("-" * 60)
print("Given Parameters:")
print(f"Top Flange Width (b1):    {b1} mm")
print(f"Bottom Flange Width (b2): {b2} mm")
print(f"Web Height (h):           {h} mm")
print(f"Flange Thickness (t_f):   {t_f} mm")
print(f"Web Thickness (t_w):      {t_w} mm")
print("-" * 60)

# 2. Calculate the Moment of Inertia (I_x) about the horizontal centroidal axis
# We use the parallel axis theorem. We'll neglect the moment of inertia of the
# flanges about their own centroidal axes as they are small for thin flanges.
# I_x for the web
I_x_web = (t_w * h**3) / 12.0
# I_x for the two flanges (Area * distance^2)
I_x_flanges = (b1 * t_f * (h / 2.0)**2) + (b2 * t_f * (h / 2.0)**2)
I_x = I_x_web + I_x_flanges

print("Calculation Steps:")
print(f"1. Moment of Inertia (I_x):")
print(f"   I_x = (t_w * h^3)/12 + t_f*b1*(h/2)^2 + t_f*b2*(h/2)^2")
print(f"   I_x = ({t_w} * {h}^3)/12 + {t_f}*{b1}*({h/2})^2 + {t_f}*{b2}*({h/2})^2")
print(f"   I_x = {I_x_web:.2f} + {I_x_flanges:.2f} = {I_x:.2f} mm^4")
print("\n2. Shear Center Offset (e) from the web centerline:")

# 3. Calculate the shear center offset 'e'
numerator = h**2 * t_f * (b1**2 + b2**2)
denominator = 8.0 * I_x
e = numerator / denominator

# 4. Print the final calculation and result
print(f"   Formula: e = (h^2 * t_f * (b1^2 + b2^2)) / (8 * I_x)")
print(f"   Substituting values:")
print(f"   e = ({h}^2 * {t_f} * ({b1}^2 + {b2}^2)) / (8 * {I_x:.2f})")
print(f"   e = ({h**2 * t_f:.2f} * ({b1**2 + b2**2:.2f})) / {denominator:.2f}")
print(f"   e = {numerator:.2f} / {denominator:.2f}")
print("-" * 60)
print(f"Result:")
print(f"The shear center is located at a distance e = {e:.2f} mm from the web centerline.")
print("This location is outside the physical cross-section, on the opposite side of the flanges.")
print("-" * 60)
