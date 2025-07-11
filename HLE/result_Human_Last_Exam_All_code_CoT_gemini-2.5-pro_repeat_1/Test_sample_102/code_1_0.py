import math

# Define the given chromium and nickel equivalent values
cr_eq = 39
ni_eq = 29

# Use the Suutala empirical formula to calculate the approximate ferrite percentage
ferrite_percentage = 79.8 * (cr_eq / ni_eq) - 68.3

# Round the calculated ferrite percentage to the nearest 10
# To round a number to the nearest 10, we divide by 10, round to the nearest integer, then multiply by 10.
rounded_ferrite_level = int(round(ferrite_percentage / 10.0) * 10)

# Print the equation with the values plugged in, showing the calculation
print(f"Using the Suutala formula to approximate ferrite level:")
print(f"Ferrite % = 79.8 * ({cr_eq} / {ni_eq}) - 68.3")
print(f"Calculated Ferrite Level = {ferrite_percentage:.2f}%")
print(f"\nRounding {ferrite_percentage:.2f} to the nearest 10 gives the final approximate level.")
print(f"Approximate Ferrite Level: {rounded_ferrite_level}")
