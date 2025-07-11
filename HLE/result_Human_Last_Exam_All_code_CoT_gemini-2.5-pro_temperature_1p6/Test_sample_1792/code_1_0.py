# Define unicode strings for the ordinal symbols for a clear output
omega = "\u03C9"
omega_1 = f"{omega}\u2081"
omega_2 = f"{omega}\u2082"
alpha_1_sym = f"\u03B1\u2081"
alpha_2_sym = f"\u03B1\u2082"
alpha_3_sym = f"\u03B1\u2083"
alpha_4_sym = f"\u03B1\u2084"
mult_dot = "\u00B7"

# The values of the coefficients in the final form, based on our derivation
alpha_1_val = omega_1
alpha_2_val = "0"
alpha_3_val = "0"
alpha_4_val = "0"

# Print the final result in the requested format.
# "Remember in the final code you still need to output each number in the final equation!"
print("The expression is expressed in the form:")
print(f"{omega_2} {mult_dot} {alpha_1_sym} + {omega_1} {mult_dot} {alpha_2_sym} + {omega} {mult_dot} {alpha_3_sym} + {alpha_4_sym}\n")

print("The determined coefficients are:")
print(f"{alpha_1_sym} = {alpha_1_val}")
print(f"{alpha_2_sym} = {alpha_2_val}")
print(f"{alpha_3_sym} = {alpha_3_val}")
print(f"{alpha_4_sym} = {alpha_4_val}\n")

print("The full equation with each term is:")
# This line prints the full equation with the values substituted.
print(f"{omega_2} {mult_dot} ({alpha_1_val}) + {omega_1} {mult_dot} ({alpha_2_val}) + {omega} {mult_dot} ({alpha_3_val}) + {alpha_4_val}\n")

print("The simplified final result is:")
print(f"{omega_2} {mult_dot} {omega_1}")