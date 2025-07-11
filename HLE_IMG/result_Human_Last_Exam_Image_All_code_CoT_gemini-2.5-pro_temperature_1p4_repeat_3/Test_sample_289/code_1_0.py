# The elevation of Descanso, CA according to GNIS is 3514 feet.
elevation_ft = 3514

# The conversion factor from feet to meters.
ft_to_m_conversion = 0.3048

# Calculate the elevation in meters.
elevation_m = elevation_ft * ft_to_m_conversion

# Print the calculation and the result.
print(f"The elevation of the unincorporated community is {elevation_ft} feet.")
print(f"To convert feet to meters, we multiply by {ft_to_m_conversion}.")
print(f"Calculation: {elevation_ft} * {ft_to_m_conversion} = {elevation_m}")
print(f"The elevation in meters is {elevation_m:.2f} meters.")