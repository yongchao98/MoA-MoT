import math

# Step 1: Define parameters of the initial container and calculate its surface area.
initial_side = 12.0
initial_surface_area = 6 * initial_side**2

# Step 2: Define the parameters for the proposed, more efficient cylindrical container.
# This design is based on a dense packing of 28 balls (4 layers of 7).
# Dimensions are chosen to be the smallest multiples of 0.5 cm that can contain the packing.

# Radius (r): A 7-ball hexagonal cluster on a 0.5cm grid requires a cylinder radius of at least 6.024 cm.
# The nearest valid dimension is 6.5 cm.
r = 6.5

# Height (h): Staggered hexagonal layers can have a vertical center-to-center separation of 3.5 cm.
# The total height for 4 layers is (3 separations * 3.5 cm) + (ball diameter of 4 cm) = 14.5 cm.
h = 14.5

# Step 3: Calculate the surface area of the new cylindrical container.
# The formula for the surface area of a cylinder is 2*pi*r^2 + 2*pi*r*h.
surface_area_base = 2 * math.pi * r**2
surface_area_side = 2 * math.pi * r * h
total_surface_area = surface_area_base + surface_area_side

# Step 4: Compare the new surface area with the original and print the result.
if total_surface_area < initial_surface_area:
    print("Yes, a more efficient container can be designed.")
    print(f"Initial Box Surface Area: 6 * {initial_side:.1f}^2 = {initial_surface_area:.2f} cm^2")
    print(f"Proposed Cylinder: r={r:.1f} cm, h={h:.1f} cm, Capacity=28 balls")
    print("\n--- Surface Area Calculation for the Cylinder ---")
    print(f"Formula: 2 * pi * r^2 + 2 * pi * r * h")
    print(f"Values: 2 * {math.pi:.4f} * {r:.1f}^2 + 2 * {math.pi:.4f} * {r:.1f} * {h:.1f}")
    print(f"Result: {surface_area_base:.2f} + {surface_area_side:.2f} = {total_surface_area:.2f} cm^2")
    
    # Format the final answer string as d[X] with the area rounded to two decimal places.
    area_val = round(total_surface_area, 2)
    description = f"cylinder r={r:.1f}, h={h:.1f}"
    final_output_string = f"{area_val:.2f}[{description}]"
    
    print("\n" + "="*20)
    print("Final Answer:")
    print(final_output_string)
    print("="*20)

else:
    # If no better solution was found, output 0.
    final_output_string = "0"
    print(final_output_string)
