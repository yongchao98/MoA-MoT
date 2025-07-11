import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The length of the rope
L = 7.0 / 2.0

# Step 1: Calculate the area in the three unobstructed quadrants (Q1, Q2, Q4).
# Each is a right triangle with legs of length L.
area_per_quadrant = 0.5 * L**2
area_q1_q2_q4 = 3 * area_per_quadrant

# Step 2: Calculate the accessible area in Quadrant 3 (x<=0, y<=0).
# We decompose the accessible region in Q3 into several disjoint geometric shapes.
# The calculations are based on finding the area of regions bounded by the house walls,
# axis-parallel lines, and the outer boundary where the rope is fully extended.

# Region D: Area in the strip -1<=x<=0, for y<=-2. This forms a trapezoid.
# Vertices: (0,-2), (-1,-2), (-1,-2.5), (0,-3.5). Height=1, Bases=1.5 and 0.5.
area_d = 0.5 * (1.5 + 0.5) * 1.0

# Region E: Area in the strip -1<=y<=0, for x<=-2. Symmetric to Region D.
area_e = 1.0

# Region F: Area in the square-like region -2<=x<=-1, y<=-2. This forms a trapezoid.
# Vertices: (-1,-2), (-2,-2), (-2,-1.5), (-1,-2.5). Height=1, Bases=0.5 and 0.5.
area_f = 0.5 * (0.5 + 0.5) * 1.0

# Region G: Area in the square-like region -2<=y<=-1, x<=-2. Symmetric to Region F.
area_g = 0.5

# Region H: Area in the square -2<=x<=-1, -2<=y<=-1.
# This 1x1 square is cut by the line x+y=-3.5, removing a small triangle.
# The cut-off triangle has area 0.5 * 0.5 * 0.5 = 0.125.
area_h = 1.0 - 0.125

# Region I: Area for x<=-2, y<=-2. This is a small triangle.
# Vertices: (-2,-2), (-2,-1.5), (-1.5,-2).
area_i = 0.125

# Total area in Q3 is the sum of these disjoint regions.
area_q3 = area_d + area_e + area_f + area_g + area_h + area_i

# Step 3: Calculate the total area.
total_area = area_q1_q2_q4 + area_q3

# Print the breakdown of the final calculation
print("The total area is the sum of the areas in the unobstructed quadrants and the accessible area in the obstructed quadrant.")
print(f"Rope Length (L) = {L}")
print(f"Area in Q1, Q2, and Q4 = 3 * (1/2 * {L}^2) = {area_q1_q2_q4}")
print(f"Area in Q3 = {area_d} + {area_e} + {area_f} + {area_g} + {area_h} + {area_i} = {area_q3}")
print(f"Total Area = {area_q1_q2_q4} + {area_q3} = {total_area}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the output to the user
print(output)