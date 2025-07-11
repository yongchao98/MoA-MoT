# The rope has a taxi-cab length R.
R = 7 / 2

# 1. Calculate the total area the horse could reach without any obstacles.
# The region defined by |x| + |y| <= R is a diamond shape.
# The formula for the area of this diamond is 2 * R^2.
total_diamond_area = 2 * R**2

# 2. Define the area of the house obstacle.
# The problem states the house is composed of three unit squares.
house_area = 3.0

# 3. Determine the final accessible area.
# In taxi-cab geometry, a path that bends around an obstacle does not necessarily
# increase in length. For this specific problem, the ability of the rope to bend
# along the house's boundary means there is no "shadow" effect that reduces the
# horse's reach.
# Therefore, the accessible region is the total diamond area minus the area
# of the house itself.
final_area = total_diamond_area - house_area

print("The horse is on a rope of taxi-cab length R.")
print(f"R = {R}")
print("\nStep 1: Calculate the maximum reachable area without obstacles.")
print("This area is a diamond shape with area = 2 * R^2.")
print(f"Maximum Area = 2 * ({R})^2 = {total_diamond_area}")
print("\nStep 2: Subtract the area of the house.")
print(f"The house area is {house_area} square units.")
print("\nStep 3: Final Calculation.")
print("The final area is the maximum area minus the house area.")
print(f"Final Area = {total_diamond_area} - {house_area} = {final_area}")