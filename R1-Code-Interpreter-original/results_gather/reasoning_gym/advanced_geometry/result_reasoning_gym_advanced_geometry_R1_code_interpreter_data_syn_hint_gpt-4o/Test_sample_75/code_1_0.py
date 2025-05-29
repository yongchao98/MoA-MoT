# Re-evaluate the slopes of the altitudes
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_AC

# Re-evaluate the intercepts of the altitudes
intercept_A = A[1] - slope_altitude_A * A[0]
intercept_B = B[1] - slope_altitude_B * B[0]

# Re-calculate the intersection of the altitudes
x_orthocenter = (intercept_B - intercept_A) / (slope_altitude_A - slope_altitude_B)
y_orthocenter = slope_altitude_A * x_orthocenter + intercept_A

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))