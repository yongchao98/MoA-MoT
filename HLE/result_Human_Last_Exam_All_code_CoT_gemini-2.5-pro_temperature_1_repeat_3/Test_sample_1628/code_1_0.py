# Define the dimensions of the asymmetric channel section
# For this calculation, we use a simplified model assuming thin walls
# and a uniform thickness (t) for all components (flanges and web).
# The thickness 't' cancels out in the final formula for the shear center location 'e'.

# b1: Width of the top flange (from web centerline to tip)
b1 = 60.0  # units: mm
# b2: Width of the bottom flange (from web centerline to tip)
b2 = 40.0  # units: mm
# h: Height of the web (distance between flange centerlines)
h = 100.0 # units: mm

# The shear center 'e' is the distance from the centerline of the web
# to the point where a transverse shear force can be applied without causing torsion.
# For a thin-walled, asymmetric channel with uniform thickness, the formula is derived
# from balancing the moment from shear flow in the flanges with the applied shear force moment.

# Formula: e = (b1^2 + b2^2) / ( (2*h)/3 + 2*(b1+b2) )

# Calculate the numerator and denominator separately for clarity
numerator = b1**2 + b2**2
denominator = (2 * h) / 3 + 2 * (b1 + b2)

# Calculate the shear center location 'e'
e = numerator / denominator

# --- Output ---
# The final answer needs to output each number in the final equation.

# Print the explanation and the formula
print("Calculation of the shear center 'e' for a thin-walled asymmetric channel section.")
print("The location 'e' is measured from the centerline of the web.")
print(f"Given dimensions: Top Flange (b1) = {b1}, Bottom Flange (b2) = {b2}, Web Height (h) = {h}\n")
print("Formula:")
print("e = (b1^2 + b2^2) / (2*h/3 + 2*(b1 + b2))\n")

# Print the calculation with substituted values
print("Calculation:")
print(f"e = ({b1}^2 + {b2}^2) / ( (2 * {h}) / 3 + 2 * ({b1} + {b2}) )")
print(f"e = ({b1**2} + {b2**2}) / ( {2*h/3:.2f} + {2*(b1+b2)} )")
print(f"e = {numerator} / {denominator:.2f}")
print(f"e = {e:.2f} mm")

print("\nThis calculation shows that the shear center is located at a specific distance")
print("from the web, outside the section, determined by the section's dimensions.")
