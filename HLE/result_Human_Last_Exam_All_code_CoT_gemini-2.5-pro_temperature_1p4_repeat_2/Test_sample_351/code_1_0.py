# Patient data
body_surface_area = 0.95  # in m^2
fluid_rate_per_bsa = 1500 # in mL/m^2/day

# Calculate total daily fluid needs
daily_fluid_needs = body_surface_area * fluid_rate_per_bsa

# Calculate hourly fluid needs
hourly_fluid_needs = daily_fluid_needs / 24

# Output the explanation and calculation
print("Calculating maintenance fluid needs using the Body Surface Area (BSA) method.")
print("The formula is: (BSA * 1500 mL/m^2) / 24 hours\n")
print("Final Equation:")
print(f"({body_surface_area} m^2 * {fluid_rate_per_bsa} mL/m^2) / 24 hours = {hourly_fluid_needs:.1f} cc/hr")
