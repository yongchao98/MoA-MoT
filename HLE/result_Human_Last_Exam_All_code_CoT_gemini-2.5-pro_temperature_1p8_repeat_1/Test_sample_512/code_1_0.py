# The dimensions of the proposed more efficient box
l = 13.5
w = 11.5
h = 10.5

# Calculate the area of each pair of faces
area_lw = l * w
area_lh = l * h
area_wh = w * h

# Calculate the total surface area of the box
surface_area = 2 * (area_lw + area_lh + area_wh)

# Create the description string for the container
description = f"box {l}x{w}x{h}"

# Print the final result in the format d[X], where d is the surface area
# and X is the container's description. The numbers from the final equation
# are used to generate the result.
print(f"{surface_area}[{description}]")