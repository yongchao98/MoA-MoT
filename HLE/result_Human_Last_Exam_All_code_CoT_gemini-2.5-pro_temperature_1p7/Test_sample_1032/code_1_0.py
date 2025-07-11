# The rope has a taxi-cab length L.
L = 7 / 2

# The area the horse can reach without any obstacles is a diamond shape
# defined by |x| + |y| <= L. The area of this diamond is 2 * L^2.
total_diamond_area = 2 * L**2

# The house consists of three unit squares, so its area is 3.
house_area = 3

# The problem can be simplified by realizing that any point reachable
# by wrapping the rope around a corner of the house is already contained
# within the main diamond-shaped area defined by |x| + |y| <= L.
#
# Therefore, the total reachable area is the area of the main diamond
# minus the area of the house, which the horse cannot enter.
reachable_area = total_diamond_area - house_area

# We will now print the equation with the calculated numbers.
# The final equation is Area = 2 * (Rope Length)^2 - (House Area)
print(f"The total reachable area is calculated by the area of the diamond minus the area of the house.")
print(f"Rope Length (L) = {L}")
print(f"Area of the diamond = 2 * L^2 = 2 * ({L})^2 = {total_diamond_area}")
print(f"Area of the house = {house_area}")
print(f"Reachable Area = {total_diamond_area} - {house_area} = {reachable_area}")
