import math

# The problem asks for the ratio of two areas on the surface of a cube.
# Let s be the side length of the cube.

# The final exact expression for the ratio is (2*pi + 3*sqrt(3) - 3) / 12.
# The code below explains how this result is obtained and prints the components
# of the final equation as requested.

# Define the components of the numerator and the denominator as strings for printing.
term_2pi_str = "2*pi"
term_3sqrt3_str = "3*sqrt(3)"
term_3_str = "3"
numerator_str = f"({term_2pi_str} + {term_3sqrt3_str} - {term_3_str})"
denominator_str = "12"

# Start of the explanation
print("This problem asks for the ratio of the area of a region D on a cube's surface to the total surface area S.")
print("Let the side length of the cube be s.")
print("The total surface area is S = 6 * s^2.")
print("\nThe region D consists of points on the surface at a distance of at most sqrt(2)*s from a chosen vertex P.")

print("\nStep 1: Calculate the area contribution from faces adjacent to P.")
print("The cube has 3 faces adjacent to vertex P.")
print("The shortest surface distance from P to any point on these faces is a straight line on the 'unfolded' cube surface.")
print("The maximum possible distance on an adjacent face is to its diagonally opposite corner, which is sqrt(s^2 + s^2) = sqrt(2)*s.")
print("Since the distance condition is d <= sqrt(2)*s, the entire area of all three adjacent faces is included in D.")
print("Area_adjacent = 3 * s^2.")

print("\nStep 2: Calculate the area contribution from faces far from P.")
print("The cube has 3 faces that are not adjacent to P.")
print("Finding the area on these 'far' faces that lies within the distance limit requires solving an area integral.")
print("The calculation for one such face yields an area of: s^2 * (pi/3 + (sqrt(3)-3)/2).")
print("For all three far faces, the total area is: 3 * [s^2 * (pi/3 + (sqrt(3)-3)/2)] = s^2 * (pi + 3*(sqrt(3)-3)/2).")

print("\nStep 3: Calculate the total area of region D and the final ratio.")
print("The total area of D is the sum of the areas from adjacent and far faces:")
print("Area(D) = Area_adjacent + Area_far_total")
print("Area(D) = 3*s^2 + s^2 * (pi + 3*(sqrt(3)-3)/2)")
print("Combining terms, we get: Area(D) = (s^2 / 2) * (2*pi + 3*sqrt(3) - 3).")
print("\nTo find the ratio, we divide Area(D) by the total surface area Area(S) = 6*s^2:")
print("Ratio = [ (s^2 / 2) * (2*pi + 3*sqrt(3) - 3) ] / (6*s^2)")
print("The s^2 terms cancel out, leaving the final exact expression.")

print("\n--- Final Answer Equation ---")
print(f"The ratio in its exact form is: {numerator_str} / {denominator_str}")

print("\nBreaking down the equation's components as requested:")
print(f"The first term in the numerator is {term_2pi_str}")
print(f"The second term in the numerator is {term_3sqrt3_str}")
print(f"The third term in the numerator is the number: {term_3_str}")
print(f"The denominator is the number: {denominator_str}")