import math

# Step 1: Define the parameters for the two metric spaces.
# X is the interval [0, L]
L = 1.0
# Y is the unit circle. We interpret "unit" to mean unit circumference,
# making it comparable to the "unit interval" [0, 1].
C = 1.0

print("This script calculates the Gromov-Hausdorff distance between the unit interval and the unit circle.")
print("-" * 20)

# Step 2: Define the spaces and their diameters.
print("Defining the metric spaces:")
print(f"Let X be the interval [0, {L}] with the standard absolute value metric.")
print(f"Let Y be the circle with circumference {C} and the intrinsic (geodesic) metric.")
print("")

diam_X = L
diam_Y = C / 2.0
print("The diameters of these spaces are:")
print(f"Diameter of X, diam(X) = {diam_X}")
print(f"Diameter of Y, diam(Y) = {diam_Y}")
print("-" * 20)

# Step 3: Calculate a lower bound for the distance.
lower_bound = 0.5 * abs(diam_X - diam_Y)
print("A lower bound for the Gromov-Hausdorff distance d_GH(X, Y) is (1/2) * |diam(X) - diam(Y)|.")
print(f"Calculation: 0.5 * |{diam_X} - {diam_Y}| = {lower_bound}")
print(f"So, d_GH(X, Y) >= {lower_bound}")
print("-" * 20)

# Step 4: Use the known formula for the exact distance.
# This also serves as an upper bound.
print("A known theorem gives the exact distance for this case.")
print("The formula for d_GH([0, L], Circle(C)) depends on whether L >= C/2.")
print(f"In our case, L = {L} and C/2 = {diam_Y}, so the condition L >= C/2 is met.")
print("The formula is: d_GH = max(C/4, L/2 - C/4)")
print("-" * 20)

# Step 5: Perform the final calculation and print the equation.
val1 = C / 4.0
val2 = L / 2.0 - C / 4.0
result = max(val1, val2)

print("Final Equation:")
print(f"d_GH = max({C}/4, {L}/2 - {C}/4)")
print(f"     = max({val1}, {L / 2.0} - {C / 4.0})")
print(f"     = max({val1}, {val2})")
print(f"     = {result}")
print("-" * 20)

print(f"The calculated distance is {result}, which matches our lower bound.")
print("Therefore, the Gromov-Hausdorff distance between the interval [0,1] and the unit circle (of circumference 1) is 1/4.")
