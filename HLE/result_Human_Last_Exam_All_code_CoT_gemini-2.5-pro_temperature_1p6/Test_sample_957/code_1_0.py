# Define the variables as strings for symbolic representation
tau_m = "\u03C4_m"  # tau_m
rho_m = "\u03C1_m"  # rho_m
i = "i"
k0 = "k\u2080" # k_0
d = "d"

# Construct the formulas from Option D
tau_numerator = f"{tau_m}^2 * e^({i}{k0}{d})"
tau_denominator = f"1 - {rho_m}^2 * e^({i}*2{k0}{d})"
tau_expression = f"({tau_numerator}) / ({tau_denominator})"

rho_numerator = f"1 - ({rho_m} - {tau_m}^2) * e^({i}*2{k0}{d}) * {rho_m}"
rho_denominator = f"1 - {rho_m}^2 * e^({i}*2{k0}{d})"
rho_expression = f"({rho_numerator}) / ({rho_denominator})"

# Print the final equations
print("The overall transmission and reflection coefficients are:")
print(f"Transmission coefficient \u03C4 = {tau_expression}")
print(f"Reflection coefficient \u03C1 = {rho_expression}")