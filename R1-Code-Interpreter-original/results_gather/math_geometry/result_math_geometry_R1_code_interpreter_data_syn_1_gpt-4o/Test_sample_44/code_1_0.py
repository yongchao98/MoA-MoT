# Given volume
volume = 36 * math.pi

# Calculate radius
r_cubed = (3 * 36) / 4
radius = r_cubed ** (1/3)

# Calculate surface area
surface_area = 4 * math.pi * (radius ** 2)

# Output the surface area in terms of pi
surface_area_in_terms_of_pi = 4 * (radius ** 2)
print(f"{surface_area_in_terms_of_pi}π")