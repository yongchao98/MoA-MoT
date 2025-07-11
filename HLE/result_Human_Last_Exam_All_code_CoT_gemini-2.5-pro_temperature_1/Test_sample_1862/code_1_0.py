# Initial parameters
t_initial = 20  # seconds
d_initial = 50  # nm
d_new = 250     # nm

# Calculate radii
r_initial = d_initial / 2
r_new = d_new / 2

# The characteristic length for lithiation time is Volume / Surface Area.
# For the initial sphere: V/A = ((4/3)*pi*r^3) / (4*pi*r^2) = r/3
ratio_initial = r_initial / 3

# For the new quarter-sphere:
# V = (1/4) * (4/3)*pi*r^3
# A = (1/4)*(4*pi*r^2) + pi*r^2 = 2*pi*r^2 (curved surface + two flat semicircles)
# V/A = ((1/3)*pi*r^3) / (2*pi*r^2) = r/6
ratio_new = r_new / 6

# Calculate the new lithiation time
# t_new = t_initial * (ratio_new / ratio_initial)
t_new = t_initial * (ratio_new / ratio_initial)

print("The relationship for the new lithiation time is:")
print("t_new = t_initial * ((V_new / A_new) / (V_initial / A_initial))")
print("Substituting the geometric ratios (V/A):")
print(f"t_new = {t_initial} s * (({r_new} nm / 6) / ({r_initial} nm / 3))")

# Calculate intermediate values for the final printed equation
calc_ratio_new = r_new / 6
calc_ratio_initial = r_initial / 3
print("Plugging in the numbers:")
print(f"t_new = {t_initial} * ({calc_ratio_new:.2f} / {calc_ratio_initial:.2f})")

# Final calculation
final_result = t_initial * (calc_ratio_new / calc_ratio_initial)
print(f"t_new = {t_initial} * {(calc_ratio_new / calc_ratio_initial):.2f}")
print(f"The final calculated lithiation time is: {final_result:.1f} seconds.")
