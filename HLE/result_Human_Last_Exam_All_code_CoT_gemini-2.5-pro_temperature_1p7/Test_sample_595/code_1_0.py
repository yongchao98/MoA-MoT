import math

# Step 1: Define triangle properties and optimal orientation
L = 18.0
# For tan(theta) = 2, sin(theta) = 2/sqrt(5), cos(theta) = 1/sqrt(5)
sqrt5 = math.sqrt(5)

# Step 2: Calculate the horizontal (dx) and vertical (dy) projections for each side
# Leg 1
dx_leg1 = L * (1/sqrt5)
dy_leg1 = L * (2/sqrt5)

# Leg 2 is perpendicular to Leg 1
dx_leg2 = L * (2/sqrt5)
dy_leg2 = L * (1/sqrt5)

# Hypotenuse vector is the sum of the two leg vectors (if oriented properly)
dx_hyp = dx_leg1 + dx_leg2
dy_hyp = dy_leg1 - dy_leg2 # using difference for one possible orientation

# Step 3: Calculate the number of squares crossed by each side
# S = 1 + floor(|dx|) + floor(|dy|)
S_leg1 = 1 + math.floor(abs(dx_leg1)) + math.floor(abs(dy_leg1))
S_leg2 = 1 + math.floor(abs(dx_leg2)) + math.floor(abs(dy_leg2))
S_hyp = 1 + math.floor(abs(dx_hyp)) + math.floor(abs(dy_hyp))

print(f"Squares crossed by Leg 1: {S_leg1}")
print(f"Squares crossed by Leg 2: {S_leg2}")
print(f"Squares crossed by Hypotenuse: {S_hyp}")

# Step 4: Calculate the total number of unique squares crossed
# k = S1 + S2 + S3 - 3 (for the 3 vertices shared between sides)
k = S_leg1 + S_leg2 + S_hyp - 3

print("\nThe total number of squares k is calculated by summing the squares crossed by each side and subtracting 3 for the overlapping squares at the vertices.")
print(f"k = {S_leg1} + {S_leg2} + {S_hyp} - 3 = {k}")
