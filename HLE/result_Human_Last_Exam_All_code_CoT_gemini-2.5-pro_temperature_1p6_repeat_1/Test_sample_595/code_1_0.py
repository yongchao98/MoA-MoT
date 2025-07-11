import math

# Define the length of the legs
leg_length = 18
hypotenuse_length = 18 * math.sqrt(2)

# Define the angle theta such that tan(theta) = 1/2
# sin(theta) = 1/sqrt(5), cos(theta) = 2/sqrt(5)
sin_theta = 1 / math.sqrt(5)
cos_theta = 2 / math.sqrt(5)

# --- Calculate crossings for the first leg ---
# Vector components are (18*cos(theta), 18*sin(theta))
leg1_dx = leg_length * cos_theta
leg1_dy = leg_length * sin_theta
leg1_crossings = math.floor(leg1_dx) + math.floor(leg1_dy)

# --- Calculate crossings for the second leg ---
# Vector components are (-18*sin(theta), 18*cos(theta))
leg2_dx = leg_length * sin_theta # abs value
leg2_dy = leg_length * cos_theta
leg2_crossings = math.floor(leg2_dx) + math.floor(leg2_dy)

# --- Calculate crossings for the hypotenuse ---
# Vector components are (18*(cos+sin), 18*(sin-cos))
hyp_dx = abs(leg_length * (cos_theta + sin_theta))
hyp_dy = abs(leg_length * (sin_theta - cos_theta))
hyp_crossings = math.floor(hyp_dx) + math.floor(hyp_dy)

# --- Total number of squares ---
k = leg1_crossings + leg2_crossings + hyp_crossings

print("Calculation for theta = arctan(1/2):")
print(f"Leg 1 projections: ({leg1_dx:.2f}, {leg1_dy:.2f})")
print(f"Leg 1 crossings: {math.floor(leg1_dx)} + {math.floor(leg1_dy)} = {leg1_crossings}")
print(f"Leg 2 projections: ({leg2_dx:.2f}, {leg2_dy:.2f})")
print(f"Leg 2 crossings: {math.floor(leg2_dx)} + {math.floor(leg2_dy)} = {leg2_crossings}")
print(f"Hypotenuse projections: ({hyp_dx:.2f}, {hyp_dy:.2f})")
print(f"Hypotenuse crossings: {math.floor(hyp_dx)} + {math.floor(hyp_dy)} = {hyp_crossings}")
print("-" * 20)
print(f"Total squares k = {leg1_crossings} + {leg2_crossings} + {hyp_crossings} = {k}")
