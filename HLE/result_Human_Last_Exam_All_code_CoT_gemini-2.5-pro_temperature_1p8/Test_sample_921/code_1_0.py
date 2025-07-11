import math

# Input values as given in the problem
angle_degrees = 40
distance_meters = 100

# Convert angle from degrees to radians for the trigonometric function
angle_radians = math.radians(angle_degrees)

# Calculate the tangent of the angle
tan_of_angle = math.tan(angle_radians)

# Calculate the height of the palace
height = distance_meters * tan_of_angle

# As requested, output each number in the final equation:
# height = distance * tan(angle_in_radians)
print(f"{height} = {distance_meters} * {tan_of_angle}")
